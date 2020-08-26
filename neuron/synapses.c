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
//! \brief Implementation of non-inlined API in synapses.h
#include "synapses.h"
#include "neuron.h"
#include "plasticity/synapse_dynamics.h"
#include <profiler.h>
#include <debug.h>
#include <spin1_api.h>
#include <utils.h>

//! if using profiler import profiler tags
#ifdef PROFILER_ENABLED
#include "profile_tags.h"
#endif //PROFILER_ENABLED

#define TIME_CONV 2147483647


//! Globals required for synapse benchmarking to work.
uint32_t  num_fixed_pre_synaptic_events = 0;

//! The number of neurons
static uint32_t n_neurons;

//! The number of synapse types
static uint32_t n_synapse_types;

//! Ring buffers to handle delays between synapses and neurons
static weight_t *ring_buffers;

//! Ring buffer size
static uint32_t ring_buffer_size;

//! Amount to left shift the ring buffer by to make it an input
static uint32_t *ring_buffer_to_input_left_shifts;

//! Count of the number of times the ring buffers have saturated
static uint32_t saturation_count = 0;

//! \brief Number of bits needed for the synapse type and index
//! \details
//! ```
//! synapse_index_bits + synapse_type_bits
//! ```
static uint32_t synapse_type_index_bits;
//! \brief Mask to pick out the synapse type and index.
//! \details
//! ```
//! synapse_index_mask | synapse_type_mask
//! ```
static uint32_t synapse_type_index_mask;
//! Number of bits in the synapse index
static uint32_t synapse_index_bits;
//! Mask to pick out the synapse index.
static uint32_t synapse_index_mask;
//! Number of bits in the synapse type
static uint32_t synapse_type_bits;
//! Mask to pick out the synapse type.
static uint32_t synapse_type_mask;


/* PRIVATE FUNCTIONS */

//#if LOG_LEVEL >= LOG_DEBUG
//! \brief get the synapse type character
//! \param[in] synapse_type: the synapse type
//! \return a single character string describing the synapse type
static inline const char *get_type_char(uint32_t synapse_type) {
    return neuron_get_synapse_type_char(synapse_type);
}
//#endif // LOG_LEVEL >= LOG_DEBUG

//! \brief Print a synaptic row.
//!
//! Only does anything when debugging.
//! \param[in] synaptic_row: The synaptic row to print
static inline void print_synaptic_row(synaptic_row_t synaptic_row) {
#if LOG_LEVEL >= LOG_DEBUG
    log_info("Synaptic row, at address %08x Num plastic words:%u\n",
            (uint32_t) synaptic_row, synapse_row_plastic_size(synaptic_row));
    if (synaptic_row == NULL) {
        return;
    }
    log_info("----------------------------------------\n");

    // Get details of fixed region
    address_t fixed_region_address = synapse_row_fixed_region(synaptic_row);
    address_t fixed_synapses =
            synapse_row_fixed_weight_controls(fixed_region_address);
    size_t n_fixed_synapses =
            synapse_row_num_fixed_synapses(fixed_region_address);
    log_info("Fixed region %u fixed synapses (%u plastic control words):\n",
            n_fixed_synapses,
            synapse_row_num_plastic_controls(fixed_region_address));

    for (uint32_t i = 0; i < n_fixed_synapses; i++) {
        uint32_t synapse = fixed_synapses[i];
        uint32_t synapse_type = synapse_row_sparse_type(
                synapse, synapse_index_bits, synapse_type_mask);

        log_info("%08x [%3d: (w: %5u (=",
                synapse, i, synapse_row_sparse_weight(synapse));
        synapses_print_weight(synapse_row_sparse_weight(synapse),
                ring_buffer_to_input_left_shifts[synapse_type]);
        log_info(
                "nA) d: %2u, %s, n = %3u)] - {%08x %08x}\n",
                synapse_row_sparse_delay(synapse, synapse_type_index_bits),
                get_type_char(synapse_type),
                synapse_row_sparse_index(synapse, synapse_index_mask),
                SYNAPSE_DELAY_MASK, synapse_type_index_bits);
    }

    // If there's a plastic region
    if (synapse_row_plastic_size(synaptic_row) > 0) {
        log_info("----------------------------------------\n");
        address_t plastic_region_address =
                synapse_row_plastic_region(synaptic_row);
        synapse_dynamics_print_plastic_synapses(
                plastic_region_address, fixed_region_address,
                ring_buffer_to_input_left_shifts);
    }

    log_info("----------------------------------------\n");
#else
    use(synaptic_row);
#endif // LOG_LEVEL >= LOG_DEBUG
}

//! \brief Print the contents of the ring buffers.
//!
//! Only does anything when debugging.
//! \param[in] time: The current timestamp
static inline void print_ring_buffers(uint32_t time) {
#if LOG_LEVEL >= LOG_DEBUG
    log_debug("Ring Buffer at %u\n", time);
    log_debug("----------------------------------------\n");
    for (uint32_t n = 0; n < n_neurons; n++) {
        for (uint32_t t = 0; t < n_synapse_types; t++) {
            // Determine if this row can be omitted
            for (uint32_t d = 0; d < (1 << SYNAPSE_DELAY_BITS); d++) {
                if (ring_buffers[synapses_get_ring_buffer_index(
                        d + time, t, n, synapse_type_index_bits,
                        synapse_index_bits)] != 0) {
                    goto doPrint;
                }
            }
            continue;
        doPrint:
            // Have to print the row
            log_debug("%3d(%s):", n, get_type_char(t));
            for (uint32_t d = 0; d < (1 << SYNAPSE_DELAY_BITS); d++) {
                log_debug(" ");
                uint32_t ring_buffer_index = synapses_get_ring_buffer_index(
                        d + time, t, n, synapse_type_index_bits,
                        synapse_index_bits);
                synapses_print_weight(ring_buffers[ring_buffer_index],
                        ring_buffer_to_input_left_shifts[t]);
            }
            log_debug("\n");
        }
    }
    log_debug("----------------------------------------\n");
#else
    use(time);
#endif // LOG_LEVEL >= LOG_DEBUG
}

//! \brief Print the neuron inputs.
//!
//! Only does anything when debugging.
static inline void print_inputs(void) {
#if LOG_LEVEL >= LOG_DEBUG
    log_debug("Inputs\n");
    neuron_print_inputs();
#endif // LOG_LEVEL >= LOG_DEBUG
}


//! \brief This is the "inner loop" of the neural simulation.
//!
//! Every spike event could cause up to 256 different weights to
//! be put into the ring buffer.
//! \param[in] fixed_region_address: The fixed region of the synaptic matrix
//! \param[in] time: The current simulation time
static inline bool process_fixed_synapses(
    address_t fixed_region_address, uint32_t payload) {

    uint32_t time  = payload;  
    //uint8_t eit =   -1;
    
    //log_info("fixed_region_address  = %u",fixed_region_address);
    // uint32_t *synaptic_words =
    //         synapse_row_fixed_weight_controls(fixed_region_address);
    // uint32_t fixed_synapse =
    //         synapse_row_num_fixed_synapses(fixed_region_address);

    // num_fixed_pre_synaptic_events += fixed_synapse;

    //log_info("fixed_synapse  = %u", fixed_synapse);


    address_t fixed_synapses =
            synapse_row_fixed_weight_controls(fixed_region_address);
    size_t n_fixed_synapses =
            synapse_row_num_fixed_synapses(fixed_region_address);
    // log_info("Fixed region %u fixed synapses (%u plastic control words):\n",
    //         n_fixed_synapses,
    //         synapse_row_num_plastic_controls(fixed_region_address));

    for (uint32_t i = 0; i < n_fixed_synapses; i++) {
        uint32_t synapse = fixed_synapses[i];
        uint32_t synapse_type = synapse_row_sparse_type(
                synapse, synapse_index_bits, synapse_type_mask);

        // log_info("%08x [%3d: (w: %5u (=",
        //         synapse, i, synapse_row_sparse_weight(synapse));
        // synapses_print_weight(synapse_row_sparse_weight(synapse),
        //         ring_buffer_to_input_left_shifts[synapse_type]);
        // log_info(
        //         "nA) d: %2u, %s, n = %3u)] - {%08x %08x}\n",
        //         synapse_row_sparse_delay(synapse, synapse_type_index_bits),
        //         get_type_char(synapse_type),
        //         synapse_row_sparse_index(synapse, synapse_index_mask),
        //         SYNAPSE_DELAY_MASK, synapse_type_index_bits);

        //log_info("synapse  = %u, synapse_type_index_bits = %u", synapse, synapse_type_index_bits);

        // Extract components from this word
        uint32_t delay =
                synapse_row_sparse_delay(synapse, synapse_type_index_bits);
        uint32_t combined_synapse_neuron_index = synapse_row_sparse_type_index(
                synapse, synapse_type_index_mask);
        uint32_t weight = synapse_row_sparse_weight(synapse);
   
        
        uint32_t neuron_index = combined_synapse_neuron_index & 255;
        
        time = payload + delay;

        //time = neuron_update_spiketime(time,neuron_index);    


        //log_info("New mc_pkt to neuron %u:  time = %u, delay = %u, weight = %f, synapse  = %u, synapse_type = %u",  neuron_index, time,delay, weight, synapse, synapse_type);

        //log_info("Time after shifting  = %u",  time);

        neuron_add_spike(time,neuron_index);

        
        uint32_t ring_buffer_index = synapses_get_ring_buffer_index_combined(
                                            time, combined_synapse_neuron_index,
                                            synapse_type_index_bits);
        //log_info("Setting ring_buffer_index  = %u for neuron_index = %u,  time = %u,  delay = %u, weight = %u", ring_buffer_index, neuron_index, time, delay, weight);
        //
        // Add weight to current ring buffer value
        //log_info("accumulation  = %u",  ring_buffers[ring_buffer_index]);
        uint32_t accumulation = ring_buffers[ring_buffer_index] + weight;
        //log_info("New accumulation  = %u",  accumulation);

        // If 17th bit is set, saturate accumulator at UINT16_MAX (0xFFFF)
        // **NOTE** 0x10000 can be expressed as an ARM literal,
        //          but 0xFFFF cannot.  Therefore, we use (0x10000 - 1)
        //          to obtain this value
        uint32_t sat_test = accumulation & 0x10000;
        if (sat_test) {
            accumulation = sat_test - 1;
            saturation_count++;
        }

        // Store saturated value back in ring-buffer
        ring_buffers[ring_buffer_index] = accumulation;        
    }


    return true;

}


//! private method for doing output debug data on the synapses
static inline void print_synapse_parameters(void) {
// only if the models are compiled in debug mode will this method contain
// said lines.
#if LOG_LEVEL >= LOG_DEBUG
    // again neuron_synapse_shaping_params has moved to implementation
    neuron_print_synapse_parameters();
#endif // LOG_LEVEL >= LOG_DEBUG
}

/* INTERFACE FUNCTIONS */
bool synapses_initialise(
        address_t synapse_params_address, uint32_t n_neurons_value,
        uint32_t n_synapse_types_value,
        uint32_t **ring_buffer_to_input_buffer_left_shifts) {
    log_debug("synapses_initialise: starting");
    n_neurons = n_neurons_value;
    n_synapse_types = n_synapse_types_value;

    // Set up ring buffer left shifts
    ring_buffer_to_input_left_shifts =
            spin1_malloc(n_synapse_types * sizeof(uint32_t));
    if (ring_buffer_to_input_left_shifts == NULL) {
        log_error("Not enough memory to allocate ring buffer");
        return false;
    }
    spin1_memcpy(
            ring_buffer_to_input_left_shifts, synapse_params_address,
            n_synapse_types * sizeof(uint32_t));
    *ring_buffer_to_input_buffer_left_shifts =
            ring_buffer_to_input_left_shifts;

    log_debug("synapses_initialise: completed successfully");
    print_synapse_parameters();

    uint32_t n_neurons_power_2 = n_neurons;
    uint32_t log_n_neurons = 1;
    if (n_neurons != 1) {
        if (!is_power_of_2(n_neurons)) {
            n_neurons_power_2 = next_power_of_2(n_neurons);
        }
        log_n_neurons = ilog_2(n_neurons_power_2);
    }

    uint32_t n_synapse_types_power_2 = n_synapse_types;
    if (!is_power_of_2(n_synapse_types)) {
        n_synapse_types_power_2 = next_power_of_2(n_synapse_types);
    }
    uint32_t log_n_synapse_types = ilog_2(n_synapse_types_power_2);

    uint32_t n_ring_buffer_bits =
            log_n_neurons + log_n_synapse_types + SYNAPSE_DELAY_BITS;
    ring_buffer_size = 1 << (n_ring_buffer_bits);

    ring_buffers = spin1_malloc(ring_buffer_size * sizeof(weight_t));
    if (ring_buffers == NULL) {
        log_error("Could not allocate %u entries for ring buffers",
                ring_buffer_size);
    }
    for (uint32_t i = 0; i < ring_buffer_size; i++) {
        ring_buffers[i] = 0;
    }

    synapse_type_index_bits = log_n_neurons + log_n_synapse_types;
    synapse_type_index_mask = (1 << synapse_type_index_bits) - 1;
    synapse_index_bits = log_n_neurons;
    synapse_index_mask = (1 << synapse_index_bits) - 1;
    synapse_type_bits = log_n_synapse_types;
    synapse_type_mask = (1 << log_n_synapse_types) - 1;
    return true;
}

uint32_t synapses_get_ring_buffer_input(uint32_t time, uint32_t neuron_index){

    //uint32_t state = spin1_irq_disable();
    // uint32_t synapse_type = synapse_row_sparse_type(
    //             synapse, synapse_index_bits, synapse_type_mask);
    uint32_t ring_buffer_index = synapses_get_ring_buffer_index(time, 0, neuron_index,
                    synapse_type_index_bits, synapse_index_bits);
    uint32_t input = synapses_convert_weight_to_input(
                            ring_buffers[ring_buffer_index],
                            ring_buffer_to_input_left_shifts[0]);
    ring_buffers[ring_buffer_index] = 0;                                        
    //log_info("Fetching ring_buffer_index = %u, input = %u",ring_buffer_index, input );
    // Re-enable the interrupts
    //spin1_mode_restore(state);

    return(input);       


}

bool synapses_process_synaptic_row(
        uint32_t payload, synaptic_row_t row) {

    // uint32_t nbPlasticElms = synapse_row_plastic_size(row);
          
    // log_info("nbPlasticElms = %u", nbPlasticElms);        


    print_synaptic_row(row);
    spin1_delay_us(1000);
    
    address_t fixed_region_address = synapse_row_fixed_region(row);

    // Process any fixed synapses
    // **NOTE** this is done after initiating DMA in an attempt
    // to hide cost of DMA behind this loop to improve the chance
    // that the DMA controller is ready to read next synaptic row afterwards
    if(process_fixed_synapses(fixed_region_address, payload)){
        return true;
    }else{
        return false;
    }

    
}

//! \brief returns the number of times the synapses have saturated their
//!        weights.
//! \return the number of times the synapses have saturated.
uint32_t synapses_get_saturation_count(void) {
    return saturation_count;
}

//! \brief returns the counters for plastic and fixed pre synaptic events
//! based on (if the model was compiled with SYNAPSE_BENCHMARK parameter) or
//! returns 0
//! \return the counter for plastic and fixed pre synaptic events or 0
uint32_t synapses_get_pre_synaptic_events(void) {
    return (num_fixed_pre_synaptic_events +
            synapse_dynamics_get_plastic_pre_synaptic_events());
}

void synapses_flush_ring_buffers(void) {
	for (uint32_t i = 0; i < ring_buffer_size; i++) {
        ring_buffers[i] = 0;
    }
}

//! \brief allows clearing of DTCM used by synapses
//! \return true if successful
bool synapses_shut_down(void) {
    sark_free(ring_buffer_to_input_left_shifts);
    sark_free(ring_buffers);
    num_fixed_pre_synaptic_events = 0;
    saturation_count = 0;
    return true;
}
