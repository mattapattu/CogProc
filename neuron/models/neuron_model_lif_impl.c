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
        return(neuron->T_refract);
    }
}

void lambda(neuron_t * neuron, key_t key, uint32_t neuron_index){
    state_t currentState  = neuron->phase;
    uint32_t nextEventTime = neuron->eot;

    if(currentState == 2){
        //clear 32nd bit if packet is spike 
        nextEventTime = nextEventTime & (~(1 << 32));
        log_info("Neuron phase = %u, nextEventTime = %u",currentState, nextEventTime );
        while (!spin1_send_mc_packet(
                        key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                    spin1_delay_us(1);
                }
    }else if(currentState == 3){
        //time  = time + neuron->tn;
        //set 32nd bit if packet is eot messg. 
        log_info("Neuron phase = %u, nextEventTime = %u",currentState, nextEventTime );
        nextEventTime = (1 << 32) | nextEventTime;
        while (!spin1_send_mc_packet(
                        key | neuron_index, nextEventTime, WITH_PAYLOAD)) {
                    spin1_delay_us(1);
                }
    }
}

//DEVS PDEVS  simulator
int32_t neuron_model_PDevs_sim(neuron_t * neuron, uint32_t threshold,  uint32_t nextSpikeTime, key_t key, uint32_t neuron_index, input_t input){
    if(neuron->tn <= neuron->eit && neuron->tn <=  nextSpikeTime ){
        //Call deltaInt()
        log_info("Process internal event  at time = %u", neuron->eit);
        neuron_model_Devs_sim(neuron, 1,nextSpikeTime, threshold, key, neuron_index, input );
        return(1);
        
    }else if(nextSpikeTime <= neuron->eit &&  nextSpikeTime < neuron->tn ){
        //Call deltaExt()
        neuron->waitCounter = 0;
        log_info("Process spike at time = %u", nextSpikeTime);
        neuron_model_Devs_sim(neuron, 2,nextSpikeTime,  threshold, key, neuron_index, input);
        neuron_model_spiketime_pop(neuron);
        return(1);
    }
    // an expected spike has been delayed
    else if(nextSpikeTime > neuron->eit){
        log_info("New event has not arrived after X clock cycles");
        if(neuron->waitCounter > 30){
            
            //log_error("On neuron = %u an expected event at eit = %u has not arrived yet",neuron_index, neuron->eit );
            log_info("New event has not arrived after waitCounter> 30. Exit simulation");
            //Add graceful exit here
        }
        neuron->waitCounter++;
        spin1_delay_us(1000);
        return(-1);
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
    if(neuron->eot == 2147483646){
        //IF next event is at time = infinity, stop PDevs while loop
        return(0);
    }else{
        return(1);
    }
}

//DEVS atomic simulator
void neuron_model_Devs_sim(neuron_t * neuron, int16_t event_type, uint32_t nextSpikeTime, uint32_t threshold, key_t key, uint32_t neuron_index, input_t input){
    //event_type 1 - Internal event
    if(event_type == 1 ){
        lambda(neuron, key, neuron_index);
        neuron->phase  = deltaInt(neuron);
        neuron->tl = neuron->tn;
        neuron->tn = neuron->tl + ta(neuron);

    }//event_type 2 - External event
    else if(event_type == 2){
        uint32_t e = nextSpikeTime  - neuron->tl;
        neuron->phase = deltaExt(neuron, nextSpikeTime, threshold, input);
        neuron->tl = nextSpikeTime;
        neuron->tn = neuron->tl + ta(neuron);
    }

}


static inline void lif_update(uint32_t time, neuron_pointer_t neuron, input_t input_this_timestep) {

    // update membrane voltage
    REAL alpha = input_this_timestep * neuron->R_membrane + neuron->V_rest;
    // update membrane voltage
    REAL V_prev = neuron_model_update_membrane_voltage(time, neuron);
    if(V_prev < neuron->V_rest){
        V_prev = neuron->V_rest;
    }
    neuron->V_membrane = alpha - (neuron->exp_TC * (alpha - V_prev));
    
    
    //log_info("Current V_membrane = %11.4k mv",  neuron->V_membrane);
    //neuron->V_membrane += input_this_timestep * neuron->R_membrane;
    log_info("New V_membrane = %11.4k mv",  neuron->V_membrane);

}

void neuron_model_eit_update(neuron_pointer_t neuron, uint32_t time){
    
    if(time < neuron->eit){
        neuron->eit = time;
    }
}

int32_t deltaExt(neuron_pointer_t neuron, uint32_t time, uint32_t threshold, input_t input) {
	//log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 2 || neuron->phase == 3){
        log_info("Ignore input as neuron is in threshold/refractory phase");
        return(neuron->phase);
    }else{
        log_info("Calling lif_update");
        lif_update(time, neuron, input);
        log_info("neuron->V_membrane = %u", neuron->V_membrane);
        if(neuron->V_membrane >= threshold){
            neuron->V_membrane = neuron->V_reset;
            return(2);
        }else{
            return(1);
        }
    }

    
}

int32_t deltaInt(neuron_pointer_t neuron) {
	
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 0 || neuron->phase == 1){
        log_info("In phase 0 or 1 = no neuron phase change ");
        return(neuron->phase);
    }else if(neuron->phase == 2){
        return(3);
    }else if(neuron->phase == 3){
        return(1);
    }

    
}

state_t neuron_model_update_membrane_voltage(uint32_t time, neuron_t *neuron) {
    
    //Check this again!!!! -> Do we update neuron membrane voltage after every state transition ( at t= tl) ?????
    uint32_t delta_t = time - neuron->tl;
    
    float exp_factor = 1;
    
    for(uint32_t k = delta_t; k > 0; k--){
 	    exp_factor = exp_factor*neuron->exp_TC;
    }
    
    if(neuron->V_membrane > neuron->V_rest) {
          neuron->V_membrane = neuron->V_membrane * (2-exp_factor); //Membrane potential is always less than 0, so decay factor > 1 : -45*(2-0.9) = -49.5 
    }
    //log_info("time = %u,last_update_time = %u, Updated V_membrane = %11.4k mv, delta_t = %u, exp_factor = %f", time, last_update_time, neuron->V_membrane, delta_t, exp_factor);
    
    return neuron->V_membrane;
}

/* void neuron_model_has_spiked(neuron_t *restrict neuron) {
    // reset membrane voltage
    neuron->V_membrane = neuron->V_reset;

    // reset refractory timer
    //neuron->refract_timer  = neuron->T_refract;
} */

state_t neuron_model_get_voltage(neuron_pointer_t neuron) {
    return neuron->V_membrane;
}

bool neuron_model_add_spike(neuron_pointer_t neuron, uint32_t  spikeTime){
   neuron->spikeCount++;
   if(neuron->spikeCount > 10){
       log_error("spiGkeCount = %u, error storing new spike at time = %u. Exiting simulation", neuron->spikeCount, spikeTime);
       for(int32_t k =0;k<10;k++){
           log_info("neuron->spike_times[k]  = %u", neuron->spike_times[k] );
       }
       return FALSE;
   }else if(spikeTime >= neuron->spike_times[neuron->spikeCount-1] ){
       log_info("Adding new spike at time = %u", spikeTime);
       neuron->spike_times[neuron->spikeCount] = spikeTime;
       return TRUE;
   }else{
       for(uint32_t i = 0; i < 9; i++){
           if(neuron->spike_times[i] < spikeTime){
               continue;
           }else{
               for(uint32_t j = 9; j >= i; j--){
                      neuron->spike_times[j] =  neuron->spike_times[j-1];      
               }
               neuron->spike_times[i]  = spikeTime;
           }
       }
       log_info("Adding new spike at time = %u", spikeTime);
       return TRUE;
   }
}


uint32_t neuron_model_spiketime_pop(neuron_pointer_t neuron){

    uint32_t nextSpike = neuron->spike_times[0];
    neuron->spike_times[0] = 0;
    for(uint32_t i = 0; i < 9; i++){
        neuron->spike_times[i] = neuron->spike_times[i+1];
    }
    neuron->spike_times[9] = 0;   
    neuron->spikeCount--;
    log_info("Popping spike at nextSpike = %u", nextSpike);
    return(nextSpike);	
}

void neuron_model_init(neuron_t *neuron){
    neuron->spikeCount = 0;
    neuron->tl = 200000;
    neuron->tn = 200000;
    neuron->eit = 200000;
    neuron->eot = 200000;
    neuron->phase = 0;
    neuron->waitCounter = 0;
    neuron->V_membrane = neuron->V_rest;
    log_info("Initializing neuron params,spikeCount = %u", neuron->spikeCount);

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
}
