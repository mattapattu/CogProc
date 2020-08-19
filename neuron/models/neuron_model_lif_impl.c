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


//! \brief simple Leaky I&F ODE
//! \param[in,out] neuron: The neuron to update
//! \param[in] V_prev: previous voltage
//! \param[in] input_this_timestep: The input to apply


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
    uint32_t nextEventTime = (uint32_t) neuron->eot;
    //log_info("lambda: neuron %u currentState = %u, nextEventTime = %u",neuron_index,  currentState, nextEventTime );
   
    if(currentState == 2){
    //clear 32nd bit if packet is spike 
    //nextEventTime = nextEventTime & (~(1 << 31));
    //log_info("Sending Spike with key = %u, neuron_index = %u, payload = %u",key,  neuron_index, nextEventTime );
    log_info("Neuron = %u has fired at time = %u",neuron_index, nextEventTime );
    neuron->lastThresholdTime = nextEventTime;
    neuron->hasSpiked  = true;
    
    }else if(currentState == 3){
        //time  = time + neuron->tn;
        //set 32nd bit if packet is eot messg. 
        
        //log_info("Sending EOT with key = %u, neuron_index = %u, time = %u",key,  neuron_index, nextEventTime );
        nextEventTime |= (1 << 31);
        
    }else if(currentState == 0||currentState == 1){
        //log_info("Sending EOT with key = %u, neuron_index = %u, time = %u",key,  neuron_index, nextEventTime );
        nextEventTime |= (1 << 31);
        //log_info("Sending EOT with key = %u, neuron_index = %u, payload = %u",key,  neuron_index, nextEventTime );
        
    }

    if(use_key && neuron->eot < INFINITY){    
        while (!spin1_send_mc_packet(
                    key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                spin1_delay_us(50);
            }
    }
    
}

int32_t neuron_model_check_next_ev(neuron_t * neuron){
    //log_info("phase = %u, tn = %u, eit = %u, spikeCount = %u", neuron->phase, neuron->tn, neuron->eit, neuron->spikeCount);
    float nextSpikeTime = neuron->spike_times[0];
    if(nextSpikeTime == INFINITY && neuron->tn == INFINITY){
        log_info("Next events at INFINITY");
        neuron->phase = 4;
        return(0);
    }else if(neuron->tn <= neuron->eit && neuron->tn <=  nextSpikeTime ){
        if(neuron->tn < INFINITY){
            log_info("tn = %f can be executed, continue PDEVS loop",neuron->tn);
            return(1);
        }else{
            return 0;
        }
        
    }else if(nextSpikeTime <= neuron->eit &&  nextSpikeTime < neuron->tn){
        if(nextSpikeTime < INFINITY){
            log_info("nextSpikeTime = %f can be executed, continue PDEVS loop",neuron->tn);
            return(1);
        }else{
            return 0;
        }
    }else{
        log_info("TN > EIT, wait for EIT update");
        log_info("tn = %f, eit = %f, nextSpikeTime = %u", neuron->tn, neuron->eit, nextSpikeTime);
        return 0;
    }
}


//DEVS PDEVS  simulator
int32_t neuron_model_PDevs_sim(neuron_t * neuron, int32_t threshold,  uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input, bool use_key){
    
    if(neuron->tn <= neuron->eit && neuron->tn <=  nextSpikeTime ){
        //Call deltaInt()
        log_info("Neuron %u PHASE %u END at tn = %f",neuron_index, neuron->phase,  neuron->tn);
        neuron_model_Devs_sim(neuron, 1,nextSpikeTime, threshold, key, neuron_index, input, use_key);
        
    }else if(nextSpikeTime <= neuron->eit &&  nextSpikeTime < neuron->tn ){
        //Call deltaExt()
        neuron->waitCounter = 0;
        log_info("Neuron %u: EXTERNAL INPUT at %u", neuron_index, nextSpikeTime);
        neuron_model_Devs_sim(neuron, 2,nextSpikeTime,  threshold, key, neuron_index, input, use_key);
        neuron->lastProcessedSpikeTime =  neuron_model_spiketime_pop(neuron);
    }
    else{
        log_info("Cannot execute any events.Check");
        log_info("tn = %f, eit = %f, nextSpikeTime = %u", neuron->tn, neuron->eit, nextSpikeTime);
        return(-1);
    }

    if(neuron->eit < neuron->tn ){
        neuron->eot = neuron->eit; 
    }else{
        neuron->eot  = neuron->tn; 
    }
    //If ext or int event has been executed,
    // check ONCE another event can be execued immediately without 
    //new pkts recvd, exit otherwise
    return(neuron_model_check_next_ev(neuron));
}

//DEVS atomic simulator
void neuron_model_Devs_sim(neuron_t * neuron, int16_t event_type, uint32_t nextSpikeTime, int32_t threshold, key_t key, uint32_t neuron_index, input_t input, bool use_key){
    //event_type 1 - Internal event
    if(event_type == 1 ){
        //log_info("Neuron %u internal event: phase %d expired at tn=%f",neuron_index, neuron->phase, neuron->tn);
        neuron->phase  = deltaInt(neuron,key,neuron_index,use_key);
        log_info("Neuron %u NEW PHASE = %u",neuron_index, neuron->phase);
        neuron->tl = neuron->tn; //UPDATE TL
        
    }//event_type 2 - External event
    else if(event_type == 2){
        neuron->phase = deltaExt(neuron, nextSpikeTime, threshold, input);
        log_info("Neuron %u NEW PHASE = %u after spike",neuron_index, neuron->phase);
        neuron->tl = (float) nextSpikeTime + deltaT;//UPDATE TL
        
    }

    if(neuron->eit != INFINITY && neuron->tl >= neuron->eit){
        log_info("Neuron %u eit = INFINITY", neuron_index);
        neuron->eit = INFINITY;
    }

    if(ta(neuron) == INFINITY){
        neuron->tn = INFINITY; //UPDATE TN
        log_info("Neuron %u TN = INFINITY", neuron_index);
    }else{
    //Next phase change = last event time + time-advance(current-phase)
        neuron->tn = neuron->tl + ta(neuron);  //UPDATE TN
        log_info("Neuron %u TN = %f", neuron_index,neuron->tn);
   }

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

bool neuron_model_eit_update(neuron_t * neuron, float time){

    if(time < INFINITY){
       neuron->eit = time;
       return true;
    }else{
        return false;
    }
    
    // if(neuron->eit == 0 || neuron->eit > time ){
    //     neuron->eit = time;
    //     //log_info("Updating neuron eit to %f",neuron->eit );
    // }
    
}

int32_t deltaExt(neuron_t * neuron, uint32_t time, int32_t threshold, input_t input) {
	//log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 2 || neuron->phase == 3){
        //log_info("Ignore input as neuron is in threshold/refractory phase");
        return(neuron->phase);
    }else{
        //log_info("external input = %f", input);
        lif_update(time, neuron, input);
        //log_info("New V_membrane after lif_update = %f, threshold = %f",  neuron->V_membrane,threshold);
        return(1);
    }
}

int32_t deltaInt(neuron_t * neuron,key_t key, uint32_t neuron_index, bool use_key ) {
	
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);
    lambda(neuron, key, neuron_index, use_key);
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
        //log_info("Neuron in phase 2, reset V_memb and change to phase 3");
        neuron->V_membrane = neuron->V_reset;
        return(3);
    }else if(neuron->phase == 3){
        return(0);
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
    neuron->eit = INFINITY;
    neuron->eot = 0;
    neuron->phase = 0;
    neuron->waitCounter = 0;
    neuron->V_membrane = neuron->V_rest;
    for(uint32_t i = 0; i < 10; i++){
        neuron->spike_times[i] = INFINITY;
    }
    //log_info("Initializing neuron params,spikeCount = %u, tl = %f", neuron->spikeCount, neuron->tl);

}

int32_t neuron_model_get_phase(neuron_t * neuron){
    return(neuron->phase);
}

void neuron_model_print_parameters(const neuron_t *neuron) {

    //log_info("V reset       = %11.4k mv", neuron->V_reset);
    //log_info("V rest        = %11.4k mv", neuron->V_rest);

    //log_info("I offset      = %11.4k nA", neuron->I_offset);
    //log_info("R membrane    = %11.4k Mohm", neuron->R_membrane);

    //log_info("exp(-ms/(RC)) = %11.4k [.]", neuron->exp_TC);

    //log_info("T refract  neuron_model_print_state_variables   = %u timesteps", neuron->T_refract);
    //log_info("V_membrane     = %f ", neuron->V_membrane);
    
}
