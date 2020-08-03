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

int32_t ta(int32_t neuron_phase){
    if(neuron_phase == 0){
        return(2147483647);
    }else if(neuron_phase == 1){
        return(2147483647);
    }else if(neuron_phase == 2){
        return(0);
    }else if(neuron_phase == 3){
        return(neuron->T_refract);
    }
}

void lambda(state_t currentState, key_t key, int32_t neuron_index, int32_t time){
    if(currentState == 2){
        while (!spin1_send_mc_packet(
                        key | neuron_index, time, WITH_PAYLOAD)) {
                    spin1_delay_us(1);
                }
    }
}

//DEVS PDEVS  simulator
void neuron_model_PDevs_sim(neuron_t * neuron, threshold_type_t *threshold_type,  int32_t nextSpikeTime, key_t key, int32_t neuron_index){
    if(neuron->tn <= neuron->eit && neuron->tn <=  nextSpikeTime ){
        //Call deltaInt()
        neuron_model_devs_sim(neuron, 1,nextSpikeTime, threshold_type, key, neuron_index );
    }else if(nextSpikeTime <= neuron->eit &&  nextSpikeTime < tn ){
        //Call deltaExt()
        neuron_model_devs_sim(neuron, 2,nextSpikeTime,  threshold_type, key, neuron_index);
    }
    int32_t lookahead = 0;

    //REDO with bitmasks
    if(neuron->eit < neuron->tn ){
        lookahead = neuron->eit; //
    }else{
        lookahead  = neuron->tn; // Infinity
    }
}

//DEVS atomic simulator
void neuron_model_Devs_sim(neuron_t * neuron, int16_t event_type, int32_t nextSpikeTime, threshold_type_t *threshold_type, key_t key){

    if(event_type == 1 ){
        lambda(neuron->phase, key, neuron_index, neuron=>tn);
        state_t neuron->phase  = deltaInt(neuron);
        neuron->tl = neuron->tn;
        neuron->tn = neuron->tl + ta(neuron->phase);

    }else if(event_type == 2){
        int32_t e = nextSpikeTime  - neuron->tl;
        state_t neuron->phase = deltaExt(neuron, e, threshold_type);
        neuron->tl = nextSpikeTime;
        neuron->tn = neuron->tl + ta(neuron->phase);
    }

}


static inline void lif_update(int32_t time, neuron_pointer_t neuron, input_t input_this_timestep) {

    // update membrane voltage
    REAL alpha = input_this_timestep * neuron->R_membrane + neuron->V_rest;
    // update membrane voltage
    REAL V_prev = neuron_model_update_membrane_voltage(time, neuron);
    neuron->V_membrane = alpha - (neuron->exp_TC * (alpha - V_prev));
    
    //log_info("Current V_membrane = %11.4k mv",  neuron->V_membrane);
    //neuron->V_membrane += input_this_timestep * neuron->R_membrane;
    log_info("New V_membrane = %11.4k mv",  neuron->V_membrane);

}

state_t deltaExt(int32_t time, 
		uint16_t num_excitatory_inputs, const input_t *exc_input,
		uint16_t num_inhibitory_inputs, const input_t *inh_input,
		input_t external_bias, neuron_t *restrict neuron, threshold_type_t *threshold_type) {
	log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 2 || neuron->phase == 3){
        log_info("Ignore input as neuron is in threshold/refractory phase")
        return(neuron->phase);
    }else{
        REAL total_exc = 0;
		REAL total_inh = 0;

		for (int i=0; i < num_excitatory_inputs; i++) {
			total_exc += exc_input[i];
		}
		for (int i=0; i< num_inhibitory_inputs; i++) {
			total_inh += inh_input[i];
		}
        // Get the input in nA
        input_t input_this_timestep =
                total_exc - total_inh + external_bias + neuron->I_offset;
	
        log_info("input_this_timestep = %12.6k, time = %u", input_this_timestep, time);

        lif_update(time, neuron, input_this_timestep);

        if(neuron->V_membrane >=threshold_type->threshold_value){
            neuron->V_membrane = neuron->V_reset;
            return(2);
        }else{
            return(1);
        }
    }

    
}

state_t deltaInt(neuron_t *restrict neuron) {
	log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    if(neuron->phase == 0 || neuron->phase == 1){
        log_info("In phase 0 or 1 = no neuron phase change ")
        return(neuron->phase);
    }else if(neuron->phase == 2){
        return(3);
    }else if(neuron->phase == 3){
        return(1);
    }

    
}

state_t neuron_model_update_membrane_voltage(int32_t time, neuron_t *neuron) {
    int32_t delta_t = time - last_update_time;
    float exp_factor = 1;
    for(int32_t k = delta_t; k > 0; k--){
 	exp_factor = exp_factor*neuron->exp_TC;
    }
    if(neuron->V_membrane > neuron->V_rest) {
          neuron->V_membrane = neuron->V_membrane * (2-exp_factor); //Membrane potential is always less than 0, so decay factor > 1 : -45*(2-0.9) = -49.5 
    }
    log_info("time = %u,last_update_time = %u, Updated V_membrane = %11.4k mv, delta_t = %u, exp_factor = %f", time, last_update_time, neuron->V_membrane, delta_t, exp_factor);
    last_update_time = time; 
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

void neuron_impl_add_spike(neuron_pointer_t neuron, int32_t  spikeTime){
   neuron->spikeCount++;
   if(neuron->spikeCount > 10){
       log_error("Error storing new spike at time = %u. Exiting simulation", spikeTime);
   }else{
       for(uint32_t i = 0; i < 10; i++){
           if(neuron->spike_times[i] < spikeTime){
               continue;
           }else{
                for(uint32_t k = 10; k > i; k--){
                    neuron->spike_times[k] = neuron->spike_times[k-1];
                }
                neuron->spike_times[i]  = spikeTime;
           }
       }
  }
   neuron->spike_times[neuron->spikeCount-1] = spikeTime;
}


int32_t neuron_model_get_next_spiketime(neuron_pointer_t neuron){
    return(neuron->spike_times[0]);	
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
