/*
 * Copyright (c) 2017-2019 The University of Manchester
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
//! \brief Inlined neuron implementation following standard component model
#ifndef _NEURON_IMPL_STANDARD_H_
#define _NEURON_IMPL_STANDARD_H_

#include "neuron_impl.h"


// Includes for model parts used in this implementation
#include <neuron/models/neuron_model.h>
#include <neuron/input_types/input_type.h>
#include <neuron/additional_inputs/additional_input.h>
#include <neuron/threshold_types/threshold_type.h>
#include <neuron/synapse_types/synapse_types.h>
#include <neuron/synapses.h>

// Further includes
#include <debug.h>

//! Indices for recording of words
enum word_recording_indices {
    //! V (somatic potential) recording index
    V_RECORDING_INDEX = 0,
    //! Gsyn_exc (excitatory synaptic conductance/current) recording index
    GSYN_EXC_RECORDING_INDEX = 1,
    //! Gsyn_inh (excitatory synaptic conductance/current) recording index
    GSYN_INH_RECORDING_INDEX = 2,
    //! Number of recorded word-sized state variables
    N_RECORDED_VARS = 3
};

//! Indices for recording of bitfields
enum bitfield_recording_indices {
    //! Spike event recording index
    SPIKE_RECORDING_BITFIELD = 0,
    //! Number of recorded bitfields
    N_BITFIELD_VARS = 1
};

// This import depends on variables defined above
#include <neuron/neuron_recording.h>

//! Array of neuron states
static neuron_t *neuron_array;

//! Input states array
static input_type_t *input_type_array;

//! Additional input array
static additional_input_t *additional_input_array;

//! Threshold states array
static threshold_type_t *threshold_type_array;

//! Global parameters for the neurons
static global_neuron_params_t *global_parameters;

//! The synapse shaping parameters
static synapse_param_t *neuron_synapse_shaping_params;

//! The number of steps to run per timestep
static uint n_steps_per_timestep;

#ifndef SOMETIMES_UNUSED
#define SOMETIMES_UNUSED __attribute__((unused))
#endif // !SOMETIMES_UNUSED

SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Initialise the particular implementation of the data
//! \param[in] n_neurons: The number of neurons
//! \return True if successful
static bool neuron_impl_initialise(uint32_t n_neurons) {
    // allocate DTCM for the global parameter details
    if (sizeof(global_neuron_params_t)) {
        global_parameters = spin1_malloc(sizeof(global_neuron_params_t));
        if (global_parameters == NULL) {
            log_error("Unable to allocate global neuron parameters"
                    "- Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for neuron array
    if (sizeof(neuron_t)) {
        neuron_array = spin1_malloc(n_neurons * sizeof(neuron_t));
        if (neuron_array == NULL) {
            log_error("Unable to allocate neuron array - Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for input type array and copy block of data
    if (sizeof(input_type_t)) {
        input_type_array = spin1_malloc(n_neurons * sizeof(input_type_t));
        if (input_type_array == NULL) {
            log_error("Unable to allocate input type array - Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for additional input array and copy block of data
    if (sizeof(additional_input_t)) {
        additional_input_array =
                spin1_malloc(n_neurons * sizeof(additional_input_t));
        if (additional_input_array == NULL) {
            log_error("Unable to allocate additional input array"
                    " - Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for threshold type array and copy block of data
    if (sizeof(threshold_type_t)) {
        threshold_type_array =
                spin1_malloc(n_neurons * sizeof(threshold_type_t));
        if (threshold_type_array == NULL) {
            log_error("Unable to allocate threshold type array - Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for synapse shaping parameters
    if (sizeof(synapse_param_t)) {
        neuron_synapse_shaping_params =
                spin1_malloc(n_neurons * sizeof(synapse_param_t));
        if (neuron_synapse_shaping_params == NULL) {
            log_error("Unable to allocate synapse parameters array"
                    " - Out of DTCM");
            return false;
        }
    }

    return true;
}

SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Add inputs to the neuron
//! \param[in] synapse_type_index: the synapse type (e.g. exc. or inh.)
//! \param[in] neuron_index: the index of the neuron
//! \param[in] weights_this_timestep: weight inputs to be added
static void neuron_impl_add_inputs(
        index_t synapse_type_index, index_t neuron_index,
        input_t weights_this_timestep) {
    // simple wrapper to synapse type input function
    synapse_param_t *parameters =
            &neuron_synapse_shaping_params[neuron_index];
    synapse_types_add_neuron_input(synapse_type_index,
            parameters, weights_this_timestep);
}

//! \brief The number of _words_ required to hold an object of given size
//! \param[in] size: The size of object
//! \return Number of words needed to hold the object (not bytes!)
static uint32_t n_words_needed(size_t size) {
    return (size + (sizeof(uint32_t) - 1)) / sizeof(uint32_t);
}

SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Load in the neuron parameters
//! \param[in] address: SDRAM block to read parameters from
//! \param[in] next: Offset of next address in store
//! \param[in] n_neurons: number of neurons
static void neuron_impl_load_neuron_parameters(
        address_t address, uint32_t next, uint32_t n_neurons) {
    log_debug("reading parameters, next is %u, n_neurons is %u ",
            next, n_neurons);

    // Read the number of steps per timestep
    n_steps_per_timestep = address[next++];
    if (n_steps_per_timestep > 1) {
        log_info("Looping over %u steps each timestep", n_steps_per_timestep);
    } else if (n_steps_per_timestep == 0) {
        log_error("bad number of steps per timestep: 0");
    }

    if (sizeof(global_neuron_params_t)) {
        log_debug("writing neuron global parameters");
        spin1_memcpy(global_parameters, &address[next],
                sizeof(global_neuron_params_t));
        next += n_words_needed(sizeof(global_neuron_params_t));
    }

    if (sizeof(neuron_t)) {
        log_debug("reading neuron local parameters");
        spin1_memcpy(neuron_array, &address[next],
                n_neurons * sizeof(neuron_t));
        next += n_words_needed(n_neurons * sizeof(neuron_t));
    }

    if (sizeof(input_type_t)) {
        log_debug("reading input type parameters");
        spin1_memcpy(input_type_array, &address[next],
                n_neurons * sizeof(input_type_t));
        next += n_words_needed(n_neurons * sizeof(input_type_t));
    }

    if (sizeof(threshold_type_t)) {
        log_debug("reading threshold type parameters");
        spin1_memcpy(threshold_type_array, &address[next],
                n_neurons * sizeof(threshold_type_t));
        next += n_words_needed(n_neurons * sizeof(threshold_type_t));
    }

    if (sizeof(synapse_param_t)) {
        log_debug("reading synapse parameters");
        spin1_memcpy(neuron_synapse_shaping_params, &address[next],
                n_neurons * sizeof(synapse_param_t));
        next += n_words_needed(n_neurons * sizeof(synapse_param_t));
    }

    if (sizeof(additional_input_t)) {
        log_debug("reading additional input type parameters");
        spin1_memcpy(additional_input_array, &address[next],
                n_neurons * sizeof(additional_input_t));
        next += n_words_needed(n_neurons * sizeof(additional_input_t));
    }

    neuron_model_set_global_neuron_params(global_parameters);

    for (index_t n = 0; n < n_neurons; n++) {
            //log_infolog_info("Initializing neuron  = %u", n);
            neuron_model_init(&neuron_array[n]);
    }

for (index_t n = 0; n < n_neurons; n++) {
        //log_info("Print neuron parameters after load");
        neuron_model_print_parameters(&neuron_array[n]);
    }
#if LOG_LEVEL >= LOG_DEBUG
    log_debug("-------------------------------------\n");
    for (index_t n = 0; n < n_neurons; n++) {
        neuron_model_print_parameters(&neuron_array[n]);
    }
    log_debug("-------------------------------------\n");
#endif // LOG_LEVEL >= LOG_DEBUG
}

// static bool neuron_impl_check_sim_end(uint32_t n_neurons){
//     bool endSim = false;
//     //bool err = false;
//     // if(!check_spiketimes_not_empty()){
//         log("in_spiketimes is empty");
//         // log("in_spiketimes is empty. Checking neuron states");
//         endSim = true;
//         for (index_t n = 0; n < n_neurons; n++) {
//             if(neuron_model_get_phase(&neuron_array[n]) == 4){
//                 endSim = endSim && true;
//             }
//         }
              
//     // }
    
//     return(endSim);
// }

static bool neuron_impl_check_sim_end(uint32_t n_neurons){
    bool endSim = true;
    bool err = false;
    if(check_spiketimes_not_empty()){
        return false;
    }    
    for (index_t n = 0; n < n_neurons; n++) {
        if(neuron_model_get_phase(&neuron_array[n]) == 4){
            endSim = endSim && true;
        }else if(neuron_model_get_phase(&neuron_array[n]) == 5){
            err = true;
            break;
        }
    }
    if(err){
        log_info("Call end sim as neuron in Error phase");
        return(err);
    }else if(endSim){
        log_info("All neurons in Idle state. Call end_sim");
        return(endSim);
    }
}


static bool neuron_impl_add_spike(index_t neuron_index, uint32_t time) {
    //log_info("Adding spike at time = %u to neuron_index  = %u", time, neuron_index);
    neuron_pointer_t neuron = &neuron_array[neuron_index];
    //log_info("neuron %u: tl = %u", neuron_index, neuron->tl);
    return(neuron_model_add_spike(neuron, time));
}

static int32_t  neuron_impl_neuron_update(uint32_t time, index_t neuron_index,
        input_t external_bias, key_t key,  bool eit, bool use_key) {
    // Get the neuron itself
    neuron_pointer_t neuron = &neuron_array[neuron_index];

    if(eit){
        float eit = (float) time;
        log_info("New eit = %f",eit );
        neuron_model_eit_update(neuron, time);
    }else{
         if(!neuron_impl_add_spike(neuron_index, time)){
            log_error("Unable to add spike to neuron %u at time = %u", neuron_index, time);
        }
        //Next input on this synapse can at the earliest occur only after curr_time + refractory_period
        //float eit = (float) time + 0.1; // Do not hardcode refraactory period, change this
        //log_info("New eit = %f",eit );
        //neuron_model_eit_update(neuron, eit); 
    }

    threshold_type_pointer_t threshold_type =
            &threshold_type_array[neuron_index];
    int32_t threshold = threshold_type->threshold_value;
    //neuron_pointer_t neuron = &neuron_array[neuron_index];

    float nextSpikeTime = neuron->spike_times[0];
    input_t input = synapses_get_ring_buffer_input(nextSpikeTime,neuron_index );
    //log_info("nextSpikeTime = %u, input = %u", nextSpikeTime, input);
    
    int32_t ret = 1;
    while(ret == 1){
        //log_info("Calling neuron_model_PDevs_sim");
        //log_info("neuron %u: tl = %u", neuron_index, neuron->tl);
        nextSpikeTime = neuron->spike_times[0];
        ret = neuron_model_PDevs_sim(neuron, threshold, nextSpikeTime, key, neuron_index, input,use_key);
        log_info("neuron_model_PDevs_sim returns %u", ret);
        // if(ret == 1){
        //     continue;
        // }else{
        //     break;
        // }
        
    }
    // if( ret == 0){
    //     //log_info("No event to process. Wait for new spike");
    //     return 0;
    // }else if(ret == -1){
    //     //log_info("New event has not arrived after X clock cycles");
    //     return -1;
    // }


    // Return the boolean to the model timestep update
    
}



SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Stores neuron parameters back into SDRAM
//! \param[out] address: the address in SDRAM to start the store
//! \param[in] next: Offset of next address in store
//! \param[in] n_neurons: number of neurons
static void neuron_impl_store_neuron_parameters(
        address_t address, uint32_t next, uint32_t n_neurons) {
    log_debug("writing parameters");

    // Skip over the steps per timestep
    next += 1;

    if (sizeof(global_neuron_params_t)) {
        log_debug("writing neuron global parameters");
        spin1_memcpy(&address[next], global_parameters,
                sizeof(global_neuron_params_t));
        next += n_words_needed(sizeof(global_neuron_params_t));
    }

    if (sizeof(neuron_t)) {
        log_debug("writing neuron local parameters");
        spin1_memcpy(&address[next], neuron_array,
                n_neurons * sizeof(neuron_t));
        next += n_words_needed(n_neurons * sizeof(neuron_t));
    }

    if (sizeof(input_type_t)) {
        log_debug("writing input type parameters");
        spin1_memcpy(&address[next], input_type_array,
                n_neurons * sizeof(input_type_t));
        next += n_words_needed(n_neurons * sizeof(input_type_t));
    }

    if (sizeof(threshold_type_t)) {
        log_debug("writing threshold type parameters");
        spin1_memcpy(&address[next], threshold_type_array,
                n_neurons * sizeof(threshold_type_t));
        next += n_words_needed(n_neurons * sizeof(threshold_type_t));
    }

    if (sizeof(synapse_param_t)) {
        log_debug("writing synapse parameters");
        spin1_memcpy(&address[next], neuron_synapse_shaping_params,
                n_neurons * sizeof(synapse_param_t));
        next += n_words_needed(n_neurons * sizeof(synapse_param_t));
    }

    if (sizeof(additional_input_t)) {
        log_debug("writing additional input type parameters");
        spin1_memcpy(&address[next], additional_input_array,
                n_neurons * sizeof(additional_input_t));
        next += n_words_needed(n_neurons * sizeof(additional_input_t));
    }
}

#if LOG_LEVEL >= LOG_DEBUG
//! \brief Print the inputs to the neurons
//! \param[in] n_neurons: The number of neurons
static inline void neuron_impl_print_inputs(uint32_t n_neurons) {
    bool empty = true;
    for (index_t i = 0; i < n_neurons; i++) {
        synapse_param_t *params = &neuron_synapse_shaping_params[i];
        empty = empty && (0 == bitsk(
                synapse_types_get_excitatory_input(params)
                - synapse_types_get_inhibitory_input(params)));
    }

    if (!empty) {
        log_info("-------------------------------------\n");

        for (index_t i = 0; i < n_neurons; i++) {
            synapse_param_t *params = &neuron_synapse_shaping_params[i];
            input_t input = synapse_types_get_excitatory_input(params)
                    - synapse_types_get_inhibitory_input(params);
            if (bitsk(input) != 0) {
                log_info("%3u: %12.6k (= ", i, input);
                synapse_types_print_input(params);
                log_info(")\n");
            }
        }
        log_info("-------------------------------------\n");
    }
}
//#if LOG_LEVEL >= LOG_DEBUG
//! \brief Print the synapse parameters of the neurons
//! \param[in] n_neurons: The number of neurons
inline void neuron_impl_print_synapse_parameters(uint32_t n_neurons) {
    log_debug("-------------------------------------\n");
    for (index_t n = 0; n < n_neurons; n++) {
        synapse_types_print_parameters(&neuron_synapse_shaping_params[n]);
    }
    log_debug("-------------------------------------\n");
}

//! \brief Get the synapse type character for a synapse type
//! \param[in] synapse_type: The synapse type
//! \return The descriptor character (sometimes two characters)
inline const char *neuron_impl_get_synapse_type_char(uint32_t synapse_type) {
    return synapse_types_get_type_char(synapse_type);
}
#endif // LOG_LEVEL >= LOG_DEBUG

#endif // _NEURON_IMPL_STANDARD_H_
