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

//! \brief simple Leaky I&F ODE
//! \param[in,out] neuron: The neuron to update
//! \param[in] V_prev: previous voltage
//! \param[in] input_this_timestep: The input to apply


void neuron_model_set_global_neuron_params(
        const global_neuron_params_t *params) {
    use(params);
    // Does Nothing - no params
}

static uint32_t ta(neuron_t * neuron){
    if(neuron->phase == 0){
        return(2147483646);
    }else if(neuron->phase == 1){
        return(2147483646);
    }else if(neuron->phase == 2){
        return(0);
    }else if(neuron->phase == 3){
        //return(neuron->T_refract);
        return(0.1); //Fix this, do not hardcode
    }
}

void lambda(neuron_t * neuron, key_t key, uint32_t neuron_index){
    state_t currentState  = neuron->phase;
    uint32_t nextEventTime = neuron->eot;

    if(currentState == 2){
        //clear 32nd bit if packet is spike 
        //nextEventTime = nextEventTime & (~(1 << 31));
        log_info("Sending Spike with key = %u, neuron_index = %u, payload = %u",key,  neuron_index, nextEventTime );
        while (!spin1_send_mc_packet(
                        key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                    spin1_delay_us(1000);
                }
    }else if(currentState == 3){
        //time  = time + neuron->tn;
        //set 32nd bit if packet is eot messg. 
        
        nextEventTime |= (1 << 31);
        log_info("Sending EOT with key = %u, neuron_index = %u, payload = %u",key,  neuron_index, nextEventTime );
        while (!spin1_send_mc_packet(
                        key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                    spin1_delay_us(1000);
                }
    }
}

int32_t neuron_model_check_pending_ev(neuron_t * neuron){
    //log_info("phase = %u, tn = %u, eit = %u, spikeCount = %u", neuron->phase, neuron->tn, neuron->eit, neuron->spikeCount);
    if(neuron->tn < 2147483646){
        log_info("Next Internal Event < Infinity, continue PDEVS loop");
        return 1;
    }else if(neuron->spikeCount > 0){
        log_info("spikeCount > 0, continue PDEVS loop");
        return 1;
    }else if(neuron->eit < 2147483646){
        log_info("earliest Input Time < Infinity, continue PDEVS loop");
        return 1;
    }else{
        log_info("no events to process, set phase to IDLE");
        neuron->phase = 4;
        return 0;
    }
}


//DEVS PDEVS  simulator
int32_t neuron_model_PDevs_sim(neuron_t * neuron, int32_t threshold,  uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input){
    if(neuron->tn <= neuron->eit && neuron->tn <=  nextSpikeTime ){
        //Call deltaInt()
        log_info("Executing internal event  at time = %u", neuron->tn);
        neuron_model_Devs_sim(neuron, 1,nextSpikeTime, threshold, key, neuron_index, input );
        
        
    }else if(nextSpikeTime <= neuron->eit &&  nextSpikeTime < neuron->tn ){
        //Call deltaExt()
        neuron->waitCounter = 0;
        log_info("Executing external event at time = %u", nextSpikeTime);
        neuron_model_Devs_sim(neuron, 2,nextSpikeTime,  threshold, key, neuron_index, input);
        neuron_model_spiketime_pop(neuron);
    }
    // an expected spike has been delayed
    else if(nextSpikeTime > neuron->eit){
        log_info("New event has not arrived after X clock cycles");
        if(neuron->waitCounter > 30){
            
            log_info("New event has not arrived after waitCounter> 30. Set neuron phase to ERR");
            //Add graceful exit here
            neuron->phase = 5;
            return(-1);
        }
        neuron->waitCounter++;
        spin1_delay_us(1000);
        
    }else{
        log_info("Unknown condition. Check");
        log_info("tn = %u, eit = %u, nextSpikeTime = %u", neuron->tn, neuron->eit, nextSpikeTime);
        return(-2);
    }
    //uint32_t lookahead = 0;

    //REDO with bitmasks
    if(neuron->eit < neuron->tn ){
        neuron->eot = neuron->eit; //
    }else{
        neuron->eot  = neuron->tn; // Infinity
    }
    // if(neuron->eot == 2147483646){
    //     //IF next event is at time = infinity, stop PDevs while loop
    //     return(0);
    // }else{
    //     return(1);
    // }
    return(neuron_model_check_pending_ev(neuron));
}

//DEVS atomic simulator
void neuron_model_Devs_sim(neuron_t * neuron, int16_t event_type, uint32_t nextSpikeTime, int32_t threshold, key_t key, uint32_t neuron_index, input_t input){
    //event_type 1 - Internal event
    if(event_type == 1 ){
        lambda(neuron, key, neuron_index);
        log_info("Internal event = phase %d expired",neuron->phase);
        neuron->phase  = deltaInt(neuron);
        log_info("New phase after deltaInt = phase %d",neuron->phase);
        
        neuron->tl = neuron->tn;
        if(ta(neuron) == 2147483646){
            neuron->tn = 2147483646;
        }else{
            //Next phase change = last event time + time-advance(current-phase)
            neuron->tn = neuron->tl + ta(neuron);
        }
        
        //log_info("Event 1 , new tl = %u, new tn = %u",neuron->tl,  neuron->tn);
        
    }//event_type 2 - External event
    else if(event_type == 2){
        uint32_t e = nextSpikeTime  - neuron->tl;
        log_info("External event = spike with neuron in phase %d",neuron->phase);
        neuron->phase = deltaExt(neuron, nextSpikeTime, threshold, input);
        log_info("New neuron phase = %d",neuron->phase);
        neuron->tl = nextSpikeTime;
        if(ta(neuron) == 2147483646){
            neuron->tn = 2147483646;
        }else{
            neuron->tn = neuron->tl + ta(neuron);
        }
        //log_info("Event 2 , new tl = %u, new tn = %u",neuron->tl,  neuron->tn);
    }

}


static inline void lif_update(uint32_t time, neuron_t * neuron, input_t input_this_timestep) {

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

void neuron_model_eit_update(neuron_t * neuron, uint32_t time){
    
    if(time < neuron->eit){
        neuron->eit = time;
    }
}

int32_t deltaExt(neuron_t * neuron, uint32_t time, int32_t threshold, input_t input) {
	//log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 2 || neuron->phase == 3){
        //log_info("Ignore input as neuron is in threshold/refractory phase");
        return(neuron->phase);
    }else{
        log_info("external input = %f", input);
        lif_update(time, neuron, input);
        log_info("New V_membrane after lif_update = %f, threshold = %f",  neuron->V_membrane,threshold);
        if(neuron->V_membrane >= -50){ //do not hard-code, change this
        //log_info("New neuron state = %d",  2);    
            return(2);
        }else{
        //log_info("New neuron state = %d",  1);        
            return(1);
        }
    }

    
}

int32_t deltaInt(neuron_t * neuron) {
	
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 0 || neuron->phase == 1){
        //log_info("Neuron in phase 0/1, no neuron phase change ");
        return(neuron->phase);
    }else if(neuron->phase == 2){
        //log_info("Neuron in phase 2, reset V_memb and change to phase 3");
        neuron->V_membrane = neuron->V_reset;
        return(3);
    }else if(neuron->phase == 3){
        return(1);
    }

    
}

int32_t neuron_model_update_membrane_voltage(uint32_t time, neuron_t *neuron) {
    
    //Check this again!!!! -> Do we update neuron membrane voltage after every state transition ( at t= tl) ?????

    uint32_t delta_t = time - neuron->tl;
    uint32_t simulation_timestep = 1000; //Redo later to read from PyNN
    uint32_t loopMax = delta_t/simulation_timestep;
    float exp_factor = neuron->exp_TC;
    
    //log_info("exp_TC  = %f, time = %u, tl = %u,  delta_t = %u, loopMax = %u, V_membrane = %f",neuron->exp_TC,time, neuron->tl,   delta_t, loopMax, neuron->V_membrane);
   
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

/* void neuron_model_has_spiked(neuron_t *restrict neuron) {
    // reset membrane voltage
    neuron->V_membrane = neuron->V_reset;

    // reset refractory timer
    //neuron->refract_timer  = neuron->T_refract;
} */

state_t neuron_model_get_voltage(neuron_t * neuron) {
    return neuron->V_membrane;
}

bool neuron_model_add_spike(neuron_t * neuron, uint32_t  spikeTime){
   neuron->spikeCount++;
   log_info("spikeCount = %u, time = %u", neuron->spikeCount, spikeTime);
   if(neuron->spikeCount > 10){
       log_error("spikeCount = %u, error storing new spike at time = %u. Exiting simulation", neuron->spikeCount, spikeTime);
       for(int32_t k =0;k<10;k++){
           log_info("neuron->spike_times[k]  = %u", neuron->spike_times[k] );
       }
       return FALSE;
   }else if(spikeTime >= neuron->spike_times[neuron->spikeCount-1] ){
       log_info("Adding new spike time = %u at the end of array", spikeTime);
       neuron->spike_times[neuron->spikeCount] = spikeTime;
       return TRUE;
   }else{
       uint32_t i;
       for(i = 0; i < 9; i++){
           if(spikeTime < neuron->spike_times[i]){
               //log_info("Adding new spike: spike_times[%u] = %u ", neuron->spike_times[i],i);
               for(uint32_t j = 9; j >= i; j--){
                      neuron->spike_times[j] =  neuron->spike_times[j-1];      
               }
               neuron->spike_times[i]  = spikeTime;
               break;
           }else{
                continue;
           }
       }
       log_info("Adding new spike at time = %u at index = %u", spikeTime,i);
       return TRUE;
   }
}


uint32_t neuron_model_spiketime_pop(neuron_t * neuron){

    uint32_t nextSpike = neuron->spike_times[0];
    for(uint32_t i = 0; i < 9; i++){
        neuron->spike_times[i] = neuron->spike_times[i+1];
    }
    neuron->spike_times[9] = 2147483646  ;   
    neuron->spikeCount--;
    log_info("Removing spike from spike_times = %u, nextspiketime = %u", nextSpike,neuron->spike_times[0]);
    return(nextSpike);	
}

void neuron_model_init(neuron_t *neuron){
    neuron->spikeCount = 0;
    neuron->tl = 0;
    neuron->tn = 2147483646;
    neuron->eit = 2147483646;
    neuron->eot = 2147483646;
    neuron->phase = 0;
    neuron->waitCounter = 0;
    neuron->V_membrane = neuron->V_rest;
    for(uint32_t i = 0; i < 10; i++){
        neuron->spike_times[i] = 2147483646;
    }
    log_info("Initializing neuron params,spikeCount = %u, tl = %u", neuron->spikeCount, neuron->tl);

}

int32_t neuron_model_get_phase(neuron_t * neuron){
    return(neuron->phase);
}

void neuron_model_print_state_variables(const neuron_t *neuron) {
    log_debug("V membrane    = %11.4k mv", neuron->V_membrane);
}

void neuron_model_print_parameters(const neuron_t *neuron) {

    log_info("V reset       = %11.4k mv", neuron->V_reset);
    log_info("V rest        = %11.4k mv", neuron->V_rest);

    log_info("I offset      = %11.4k nA", neuron->I_offset);
    log_info("R membrane    = %11.4k Mohm", neuron->R_membrane);

    log_info("exp(-ms/(RC)) = %11.4k [.]", neuron->exp_TC);

    log_info("T refract     = %u timesteps", neuron->T_refract);
    log_info("V_membrane     = %f ", neuron->V_membrane);
    
}
