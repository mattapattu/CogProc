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

static int32_t last_update_time = 0;

//! \brief simple Leaky I&F ODE
//! \param[in,out] neuron: The neuron to update
//! \param[in] V_prev: previous voltage
//! \param[in] input_this_timestep: The input to apply

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


void neuron_model_set_global_neuron_params(
        const global_neuron_params_t *params) {
    use(params);
    // Does Nothing - no params
}

state_t neuron_model_state_update(int32_t time,
		uint16_t num_excitatory_inputs, const input_t *exc_input,
		uint16_t num_inhibitory_inputs, const input_t *inh_input,
		input_t external_bias, neuron_t *restrict neuron) {
	log_info("Exc 1: %12.6k", exc_input[0]);
	//log_info("Inh 1: %12.6k, Inh 2: %12.6k", inh_input[0], inh_input[1]);

    // If outside of the refractory period
    if (neuron->refract_timer <= 0) {
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

    } else {
        // countdown refractory timer
        neuron->refract_timer--;
    }
    return neuron->V_membrane;
}


void neuron_model_has_spiked(neuron_t *restrict neuron) {
    // reset membrane voltage
    neuron->V_membrane = neuron->V_reset;

    // reset refractory timer
    neuron->refract_timer  = neuron->T_refract;
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

state_t neuron_model_get_voltage(neuron_pointer_t neuron) {
    return neuron->V_membrane;
}

void neuron_model_set_spike_time(neuron_pointer_t neuron, int32_t  spikeTime){
   neuron->spikeCount++;
   if(neuron->spikeCount > 10){
       log_error("Error storing new spike at time = %u. Exiting simulation", spike_time);
   }else{
       for(uint32_t i = 0; i < 10; i++){
           if(neuron->spike_time[i] < spikeTime){
               continue;
           }else{
                for(uint32_t k = 10; k > i; k--){
                    neuron->spike_time[k] = neuron->spike_time[k-1];
                }
                neuron->spike_time[i]  = spikeTime;
           }
       }
  }
   neuron->spike_time[neuron->spikeCount-1] = spikeTime;
}


int32_t neuron_model_get_next_spiketime(neuron_pointer_t neuron){
    return(neuron->spike_time[0]);	
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
