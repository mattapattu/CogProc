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

#include <debug.h>

#define INFINITY  32767
#define deltaT 1 

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
    float timeadvance = 0;

    if(neuron->phase == 0 || neuron->phase == 1 ){
        timeadvance = INFINITY;
    }
    // else if(neuron->phase == 1){
    //     timeadvance = INFINITY;
    // }
    else if(neuron->phase == 2){
        timeadvance = 0;
    }else if(neuron->phase == 3){
        //return(neuron->T_refract);
        timeadvance = 0.1;  //Fix this, do not hardcode
    }
    return(timeadvance);

}

static float neuron_model_update_membrane_voltage(float time, neuron_t *neuron) {
    
    //Check this again!!!! -> Do we update neuron membrane voltage after every state transition ( at t= tl) ?????

    uint32_t delta_t = time - neuron->tl;
    //uint32_t simulation_timestep = 1000; //Redo later to read from PyNN
    //uint32_t loopMax = delta_t/simulation_timestep;
    float exp_factor = neuron->exp_TC;

    if(neuron->V_membrane > neuron->V_rest) {
          for(uint32_t k = delta_t; k > 1; k--){
 	            exp_factor = exp_factor*neuron->exp_TC;
           }  
          //log_info("exp_factor = %f, V_membrane = %f", exp_factor, neuron->V_membrane); 
          neuron->V_membrane = neuron->V_membrane * (2-exp_factor); //Membrane potential is always less than 0, so decay factor > 1 : -45*(2-0.9) = -49.5 
          //log_info("Updated V_membrane = %f, delta_t = %u, exp_factor = %f", neuron->V_membrane, delta_t, exp_factor);
    }
   
    return(neuron->V_membrane);
}

static void lambda(neuron_t * neuron, key_t key, uint32_t neuron_index, bool use_key){
    //uint16_t currentState  = neuron->phase;
    uint32_t nextEventTime = (uint32_t) neuron->eot;
    log_info("lambda: neuron %u phase = %u, nextEventTime = %u",neuron_index,  neuron->phase, nextEventTime );
    if(use_key){
        if(neuron->phase == 2){
        //clear 32nd bit if packet is spike 
        //nextEventTime = nextEventTime & (~(1 << 31));
         log_info("Sending Spike with key = %u, neuron_index = %u, payload = %u",key,  neuron_index, nextEventTime );
        }
        else if(neuron->phase == 0||neuron->phase == 1 || neuron->phase == 3){
          
            nextEventTime |= (1 << 31);
            log_info("Phase = %u, Sending EOT with key = %u, neuron_index = %u, payload = %u",neuron->phase, key,  neuron_index, nextEventTime );
           
        }
        while (!spin1_send_mc_packet(
                        key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                    spin1_delay_us(1000);
        }
    }
    
}



static void lif_update(float time, neuron_t * neuron, input_t input_this_timestep) {

    // update membrane voltage
    REAL alpha = input_this_timestep * neuron->R_membrane + neuron->V_rest;
    // update membrane voltage
    REAL V_prev = neuron_model_update_membrane_voltage(time, neuron);
    neuron->V_membrane = alpha - (neuron->exp_TC * (alpha - V_prev));
    
    

}

static int32_t deltaExt(neuron_t * neuron, uint32_t time, int32_t threshold, input_t input) {

    if(neuron->phase == 2 || neuron->phase == 3){
        //log_info("Ignore input as neuron is in threshold/refractory phase");
        return(neuron->phase);
    }else{
        //log_info("external input = %f", input);
        lif_update(time, neuron, input);
        //log_info("New V_membrane after lif_update = %f, threshold = %f",  neuron->V_membrane,threshold);
        if(neuron->V_membrane >= -50){ //do not hard-code, change this
        //log_info("New neuron state = %d",  2);    
            return(2);
        }else{
        //log_info("New neuron state = %d",  1);        
            return(1);
        }
    }

    
}

static int32_t deltaInt(neuron_t * neuron) {
	
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

static uint32_t neuron_model_spiketime_pop(neuron_t * neuron){

    uint32_t nextSpike = neuron->spike_times[0];
    for(uint8_t i = 0; i < 5; i++){
        neuron->spike_times[i] = neuron->spike_times[i+1];
    }
    neuron->spike_times[4] = INFINITY;   
    neuron->spikeCount--;
    //log_info("Removing spike from spike_times = %f, spikeCount = %u", nextSpike,neuron->spikeCount);
    return(nextSpike);	
}

//DEVS atomic simulator
static void neuron_model_Devs_sim(neuron_t * neuron, int16_t event_type, uint32_t nextSpikeTime, int32_t threshold, key_t key, uint32_t neuron_index, input_t input, bool use_key){
    //event_type 1 - Internal event
    if(event_type == 1 ){
        lambda(neuron, key, neuron_index, use_key);
        log_info("Neuron %u internal event: phase %d expired at tn=%f",neuron_index, neuron->phase, neuron->tn);
        neuron->phase  = deltaInt(neuron);
        log_info("New phase after deltaInt = phase %d",neuron->phase);
        
        neuron->tl = neuron->tn;
        if(ta(neuron) >= INFINITY){
            neuron->tn = INFINITY;
            log_info("Neuron %u: update tn = INFINITY", neuron_index);
        }else{
            //Next phase change = last event time + time-advance(current-phase)
            neuron->tn = neuron->tl + ta(neuron);
            log_info("Neuron %u: update tn = %f", neuron_index,neuron->tn);
        }
        
        log_info("Event 1 , new tl = %u, new tn = %u",neuron->tl,  neuron->tn);
        
    }//event_type 2 - External event
    else if(event_type == 2){
        //uint32_t e = nextSpikeTime  - neuron->tl;
        log_info("Neuron %u external event: at time = %f in phase %d",neuron_index, nextSpikeTime, neuron->phase);
        neuron->phase = deltaExt(neuron, nextSpikeTime, threshold, input);
        log_info("New neuron phase = %d",neuron->phase);
        neuron->tl = (float) nextSpikeTime + deltaT; //deltaExt lasts for deltaT milliseconds.
        if(ta(neuron) >= INFINITY){
            neuron->tn = INFINITY;
            log_info("Neuron %u: update tn = INFINITY", neuron_index);
        }else{
            neuron->tn = neuron->tl + ta(neuron);
            log_info("Neuron %u: update tn = %f", neuron_index,neuron->tn);
        }
        log_info("Event 2 , new tl = %u, new tn = %u",neuron->tl,  neuron->tn);
    }

}

int32_t neuron_model_check_pending_ev(neuron_t * neuron){
    //log_info("phase = %u, tn = %u, eit = %u, spikeCount = %u", neuron->phase, neuron->tn, neuron->eit, neuron->spikeCount);
    if(neuron->tn < INFINITY || neuron->spike_times[0] < INFINITY ){
        //log_info("tn = %f < INFINITY, repeat PDEVS",neuron->tn);
        log_info("tn = %f, spikeCount = %u, repeat PDEVS", neuron->tn, neuron->spikeCount);
        return 1;
    }
    // else if(neuron->spike_times[0] < INFINITY){
    //     log_info("spikeCount = %u > 0, repeat PDEVS", neuron->spikeCount);
    //     return 1;
    // }
    else{
        log_info("set phase to IDLE");
        neuron->phase = 4;
        return 0;
    }
}


//DEVS PDEVS  simulator
int32_t neuron_model_PDevs_sim(neuron_t * neuron, int32_t threshold,  uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input, bool use_key){
    if(neuron->tn <= neuron->eit && neuron->tn <=  nextSpikeTime ){
        //Call deltaInt()
        log_info("Executing tn =  phase %u expired at tn = %u",neuron->phase,  neuron->tn);
        neuron_model_Devs_sim(neuron, 1,nextSpikeTime, threshold, key, neuron_index, input, use_key);
        
        
    }else if(nextSpikeTime <= neuron->eit &&  nextSpikeTime < neuron->tn ){
        //Call deltaExt()
        log_info("Executing ext. event at time = %u", nextSpikeTime);
        neuron_model_Devs_sim(neuron, 2,nextSpikeTime,  threshold, key, neuron_index, input, use_key);
        neuron->lastProcessedSpikeTime = neuron_model_spiketime_pop(neuron);
    }else if(nextSpikeTime == INFINITY && neuron->tn == INFINITY){
        log_info("Stop nextSpikeTime & tn = INFINITY");
    }
    else{
        //log_info("Cannot process nextspike or the next internal event");
        log_info("Unknown cond: tn = %f, eit = %f, nextSpikeTime = %f, spikeCount = %u", neuron->tn, neuron->eit, nextSpikeTime, neuron->spikeCount);
        return(-1);
    }

    if(neuron->eit < neuron->tn ){
        neuron->eot = neuron->eit; 
    }else{
        neuron->eot  = neuron->tn; 
    }
    return(neuron_model_check_pending_ev(neuron));
}



void neuron_model_eit_update(neuron_t * neuron, float time){

    if(neuron->eit == 0 || neuron->eit > time ){
        neuron->eit = time;
        //log_info("Updating neuron eit to %f",neuron->eit );
    }
    
}


bool neuron_model_add_spike(neuron_t * neuron, uint32_t  spikeTime){
   neuron->spikeCount++;
   if(neuron->spikeCount > 5){
       log_error("spikeCount = %u, Exiting simulation", neuron->spikeCount);
       return FALSE;
   }else if(spikeTime >= neuron->spike_times[neuron->spikeCount-1] ){
       log_info("New spike: time = %f , index = 4, spikeCount = %u", spikeTime,neuron->spikeCount);
       neuron->spike_times[neuron->spikeCount] = spikeTime;
       return TRUE;
   }else{
       uint8_t i;
       for(i = 0; i < 4; i++){
           if(spikeTime < neuron->spike_times[i]){
               for(uint8_t j = 4; j >= i; j--){
                      neuron->spike_times[j] =  neuron->spike_times[j-1];      
               }
               neuron->spike_times[i]  = spikeTime;
               break;
           }else{
                continue;
           }
       }
       log_info("New spike: time = %f, index = %u, spikeCount = %u", spikeTime,i,neuron->spikeCount);
       return TRUE;
   }
}




void neuron_model_init(neuron_t *neuron){
    neuron->spikeCount = 0;
    neuron->tl = 0;
    neuron->tn = INFINITY;
    neuron->eit = INFINITY;
    neuron->eot = 0;
    neuron->phase = 0;
    neuron->V_membrane = neuron->V_rest;
    neuron->lastProcessedSpikeTime = 0;
    for(uint32_t i = 0; i < 5; i++){
        neuron->spike_times[i] = INFINITY;
    }
    //log_info("Initializing neuron params,spikeCount = %u, tl = %f", neuron->spikeCount, neuron->tl);

}

uint8_t neuron_model_get_phase(neuron_t * neuron){
    return(neuron->phase);
}

/* void neuron_model_print_parameters(const neuron_t *neuron) {

     log_info("V reset       = %11.4k mv", neuron->V_reset);
    log_info("V rest        = %11.4k mv", neuron->V_rest);

    log_info("I offset      = %11.4k nA", neuron->I_offset);
    log_info("R membrane    = %11.4k Mohm", neuron->R_membrane);

    log_info("exp(-ms/(RC)) = %11.4k [.]", neuron->exp_TC);

    log_info("T refract  neuron_model_print_state_variables   = %u timesteps", neuron->T_refract);
    log_info("V_membrane     = %f ", neuron->V_membrane);
    
} */
