 /*Copyright (c) 2017-2019 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! \file
//! \brief Leaky Integrate and Fire neuron implementation
#include "neuron_model_lif_impl.h"
#include <neuron/threshold_types/threshold_type.h>


#include <debug.h>

//static int32_t last_update_time = 0;
#define INFINITY  2147483646
#define deltaT 1 // deltaT = 1000ms
#define DELAY 14


//! \brief simple Leaky I&F ODE
//! \param[in,out] neuron: The neuron to update
//! \param[in] V_prev: previous voltage
//! \param[in] input_this_timestep: The input to apply

static uint32_t simulation_ticks;

void neuron_model_set_global_neuron_params(
        const global_neuron_params_t *params) {
    use(params);
    // Does Nothing - no params
}

static float ta(neuron_t * neuron){
    if(neuron->phase == 0){
        return(INFINITY);
    }else if(neuron->phase == 1 && neuron->V_membrane < -50){
        return(INFINITY);
    }else if(neuron->phase == 1 && neuron->V_membrane >= -50){
        return(0);
    }
    else if(neuron->phase == 2){
        return(0);
    }else if(neuron->phase == 3){
        //return(neuron->T_refract);
        return(0.1); //Fix this, do not hardcode
    }
    else{
        log_info("Unknown Neuron PHASE = %u. Check",  neuron->phase);
         return(-1);
    }
}

static float neuron_model_update_membrane_voltage(float time, neuron_t *neuron) {
    
    //Check this again!!!! -> Do we update neuron membrane voltage after every state transition ( at t= tl) ?????

    float delta_t = time - neuron->tl;
    uint32_t simulation_timestep = 1000; //Redo later to read from PyNN
    uint32_t loopMax = delta_t/simulation_timestep;
    float exp_factor = neuron->exp_TC;

    if(neuron->V_membrane > neuron->V_rest) {
          for(uint32_t k = loopMax; k > 1; k--){
 	            exp_factor = exp_factor*neuron->exp_TC;
           }  
          //log_info("exp_factor = %f, V_membrane = %f", exp_factor, neuron->V_membrane); 
          neuron->V_membrane = neuron->V_membrane * (2-exp_factor); //Membrane potential is always less than 0, so decay factor > 1 : -45*(2-0.9) = -49.5 
          //log_info("Updated V_membrane = %f, delta_t = %u, exp_factor = %f", neuron->V_membrane, delta_t, exp_factor);
    }
   
    return(neuron->V_membrane);
}

static void lambda(neuron_t * neuron, key_t key, uint32_t neuron_index, bool use_key){
    uint16_t currentState  = neuron->phase;
    uint32_t nextEventTime = (uint32_t) neuron->tn;
    //log_info("lambda: neuron %u currentState = %u, nextEventTime = %u",neuron_index,  currentState, nextEventTime );
   
    if(currentState == 2){
        neuron->lastThresholdTime = nextEventTime;
        neuron->hasSpiked  = true;
    }

    if(use_key && nextEventTime < INFINITY){    
        log_info("lambda: Neuron %u  sending spike at t = %u",neuron_index,  nextEventTime );
        while (!spin1_send_mc_packet(
                    key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                spin1_delay_us(1);
            }
    }
    
}


//DEVS PDEVS  simulator
int32_t neuron_model_PDevs_sim(neuron_t * neuron, int32_t threshold,  uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input, bool use_key){
    
     
    //Calling Basic DEVS
    if(neuron->tn <=  nextSpikeTime ){
        //Call deltaInt()
        log_info("Neuron %u PHASE %u END at tn = %f",neuron_index, neuron->phase,  neuron->tn);
        neuron_model_Devs_sim(neuron, 1,nextSpikeTime, key, neuron_index, input, use_key);
        
    }else if( nextSpikeTime < neuron->tn ){
        //Call deltaExt()
        neuron->waitCounter = 0;
        log_info("Neuron %u: EXTERNAL INPUT at %u", neuron_index, nextSpikeTime);
        if(neuron->tl <= nextSpikeTime|| neuron->tl == 0){
            neuron_model_Devs_sim(neuron, 2,nextSpikeTime,  key, neuron_index, input, use_key);
        }else {
            log_error("Causality Error: nextSpikeTime = %u, TL = %f in  neuron %u", nextSpikeTime, neuron->tl, neuron_index);
        }
        
        neuron->lastProcessedSpikeTime =  neuron_model_spiketime_pop(neuron);
    }else{
        log_info("Cannot execute any events.Check");
        log_info("Check: tn = %f, nextSpikeTime = %u, tl = %f", neuron->tn, nextSpikeTime, neuron->tn);
        return(-1);
    }

    //If ext or int event has been executed,
    // check ONCE another event can be execued immediately without 
    //new pkts recvd, exit otherwise

    if(nextSpikeTime >= INFINITY && neuron->tn >= INFINITY){
        //log_info("Next events at INFINITY");
        return(0);
    }else{
        return(1);
    }

    //return(neuron_model_check_next_ev(neuron));
}

//DEVS atomic simulator
void neuron_model_Devs_sim(neuron_t * neuron, int16_t event_type, uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input, bool use_key){
    //event_type 1 - Internal event
    if(event_type == 1 ){
        //log_info("Neuron %u internal event: phase %d expired at tn=%f",neuron_index, neuron->phase, neuron->tn);
        if(neuron->tn < INFINITY){
            neuron->tl = neuron->tn;
        }
        
        neuron->phase  = deltaInt(neuron,key,neuron_index,use_key);
         //UPDATE TL
        //log_info("Neuron %u NEW PHASE = %u, TL = %f",neuron_index, neuron->phase, neuron->tl);
        
    }//event_type 2 - External event
    else if(event_type == 2){
        neuron->phase = deltaExt(neuron, nextSpikeTime, input);
        neuron->tl = (float) nextSpikeTime + deltaT;//UPDATE TL
        log_info("Neuron %u NEW PHASE = %u, TL = %f after spike",neuron_index, neuron->phase, neuron->tl);
    }


    
    if(ta(neuron) == INFINITY){
        neuron->tn = INFINITY; //UPDATE TN
        //log_info("Neuron %u TN = INFINITY", neuron_index);
    }else{
    //Next phase change = last event time + time-advance(current-phase)
        neuron->tn = neuron->tl + ta(neuron);  //UPDATE TN
   }

}

void neuron_set_simulation_ticks(uint32_t time){
    log_info("simulation_ticks = %u", time);
    simulation_ticks = time;
}


static inline void lif_update(float time, neuron_t * neuron, input_t input_this_timestep) {

    // update membrane voltage
    REAL alpha = input_this_timestep * neuron->R_membrane + neuron->V_rest;
    //log_info("alpha = %u, time = %u, tl = %u",  alpha, time, neuron->tl);
    // update membrane voltage
    REAL V_prev = neuron_model_update_membrane_voltage(time, neuron);
    //log_info("V_prev = %f",  V_prev);
    // if(V_prev < neuron->V_rest){
    //     V_prev = neuron->V_rest;
    // }
    neuron->V_membrane = alpha - (neuron->exp_TC * (alpha - V_prev));
    
    

}


int32_t deltaExt(neuron_t * neuron, uint32_t time, input_t input) {
	//log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 2 || neuron->phase == 3){
        //log_info("Ignore input as neuron is in threshold/refractory phase");
        return(neuron->phase);
    }else if(neuron->phase == 0 || neuron->phase == 1){
        //log_info("external input = %f", input);
        lif_update(time, neuron, input);
        return(1);
    }else{
        log_info("Unknown Neuron PHASE = %u. Check",  neuron->phase);
        return(0);
    }

    
}

uint16_t deltaInt(neuron_t * neuron,key_t key, uint32_t neuron_index, bool use_key ) {
	
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);
    
    
    if(neuron->phase == 0 ){
        //log_info("Neuron in phase 0/1, no neuron phase change ");
        return(neuron->phase);
    }else if(neuron->phase == 1){

        if(neuron->V_membrane >= -50){ //do not hard-code, change this
        //log_info("New neuron state = %d",  2);    
            return(2);
        }else{
        //log_info("New neuron state = %d",  1);        
            return(1);
        }
    }else if(neuron->phase == 2){
        lambda(neuron, key, neuron_index, use_key);
        log_info("Neuron in phase 2, reset V_memb and change to phase 3");
        neuron->V_membrane = neuron->V_reset;
        return(3);
    }else if(neuron->phase == 3){

        log_info("Neuron %u TL = %f. Setting to phase 0", neuron_index, neuron->tl);
        return(0);
    }
    else{
        log_info("Unknown Neuron %u PHASE = %u. Check", neuron_index, neuron->phase);
        return(5);
    }

    
}

bool neuron_model_add_spike(neuron_t * neuron, uint32_t  spikeTime){
   neuron->spikeCount++;
   if(neuron->spikeCount > 10){
       log_error("spikeCount = %u, error storing new spike at time = %u. Exiting simulation", neuron->spikeCount, spikeTime);
       return FALSE;
   }else if(spikeTime >= neuron->spike_times[neuron->spikeCount-1] ){
       //log_info("Adding new spike time = %f at the end of array", spikeTime);
       neuron->spike_times[neuron->spikeCount] = spikeTime;
       return TRUE;
   }else{
       uint32_t i;
       for(i = 0; i < 9; i++){
           if(spikeTime < neuron->spike_times[i]){
               for(uint32_t j = 9; j >= i; j--){
                      neuron->spike_times[j] =  neuron->spike_times[j-1];      
               }
               neuron->spike_times[i]  = spikeTime;
               break;
           }else{
                continue;
           }
       }
       //log_info("Adding new spike at time = %f at index = %u", spikeTime,i);
       return TRUE;
   }
}


uint32_t neuron_model_spiketime_pop(neuron_t * neuron){

    uint32_t nextSpike = neuron->spike_times[0];
    for(uint32_t i = 0; i < 9; i++){
        neuron->spike_times[i] = neuron->spike_times[i+1];
    }
    neuron->spike_times[9] = INFINITY  ;   
    neuron->spikeCount--;
    //log_info("Removing spike from spike_times = %f, nextspiketime = %f", nextSpike,neuron->spike_times[0]);
    return(nextSpike);	
}

void neuron_model_init(neuron_t *neuron){
    neuron->spikeCount = 0;
    neuron->tl = 0;
    neuron->tn = INFINITY;
    neuron->phase = 0;
    neuron->waitCounter = 0;
    neuron->V_membrane = neuron->V_rest;
    for(uint32_t i = 0; i < 10; i++){
        neuron->spike_times[i] = INFINITY;
    }
    //log_info("Initializing neuron params,spikeCount = %u, tl = %f", neuron->spikeCount, neuron->tl);

}

uint16_t neuron_model_get_phase(neuron_t * neuron){
    return(neuron->phase);
}

// bool neuron_model_check(neuron_t * neuron){
//     if(neuron->phase == 4){
//         return(true);
//     }else{
//         return(false);
//     }
    
// }

bool neuron_model_check_sim_continue(neuron_pointer_t neuron){

    if(neuron->tl < simulation_ticks){
        log_info("neuron TL = %f < simulation_ticks, continue",neuron->tl);
        return true;
    }else{
        log_info("neuron TL = %f > simulation_ticks, exit",neuron->tl);
        return false;
    }
}

void neuron_model_print_parameters(const neuron_t *neuron) {

    log_debug("V reset       = %11.4k mv", neuron->V_reset);
    log_debug("V rest        = %11.4k mv", neuron->V_rest);

    log_debug("I offset      = %11.4k nA", neuron->I_offset);
    log_debug("R membrane    = %11.4k Mohm", neuron->R_membrane);

    log_debug("exp(-ms/(RC)) = %11.4k [.]", neuron->exp_TC);

    log_debug("T refract  neuron_model_print_state_variables   = %u timesteps", neuron->T_refract);
    log_debug("V_membrane     = %f ", neuron->V_membrane);
    
}
