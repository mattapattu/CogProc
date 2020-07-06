// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/delay_extension/delay_extension.c

/*
 * Copyright (c) 2017-2019 The University of Manchester
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! \file
//! \brief Implementation of delay extensions

#include "delay_extension.h"

#include <common/neuron-typedefs.h>
#include <common/in_spikes.h>
#include <bit_field.h>
#include <data_specification.h>
#include <debug.h>
#include <simulation.h>
#include <spin1_api.h>

//! values for the priority for each callback
enum delay_extension_callback_priorities {
    MC_PACKET = -1, //!< multicast packet reception uses FIQ
    SDP = 0,        //!< SDP handling is highest ordinary priority
    USER = 1,       //!< User interrupt is next (clearing packet received queue)
    DMA = 2,        //!< DMA complete handling is next
    TIMER = 3,      //!< Regular timer tick handling is lowest priority
};

//! Structure of the provenance data
struct delay_extension_provenance {
    //! Number of input spikes
    uint32_t n_packets_received;
    //! Number of spikes transferred via queue
    uint32_t n_packets_processed;
    //! Number of spikes added to delay processing
    uint32_t n_packets_added;
    //! Number of spikes sent
    uint32_t n_packets_sent;
    //! Number of circular buffer overflows (spikes internally dropped)
    uint32_t n_buffer_overflows;
    //! Number of times we had to back off because the comms hardware was busy
    uint32_t n_delays;
};

// Globals
//! Base multicast key for sending messages
static uint32_t key = 0;
//! Key for receiving messages
static uint32_t incoming_key = 0;
//! Mask for ::incoming_key to say which messages are for this program
static uint32_t incoming_mask = 0;
//! \brief Mask for key (that matches ::incoming_key/::incoming_mask) to extract
//! the neuron ID from it
static uint32_t incoming_neuron_mask = 0;

//! Number of neurons supported.
static uint32_t num_neurons = 0;

//! Simulation time
static uint32_t time = UINT32_MAX;
//! Simulation speed
static uint32_t simulation_ticks = 0;
//! True if we're running forever
static uint32_t infinite_run;

//! \brief The spike counters, as a 2D array
//! ```
//! spike_counters[time_slot][neuron_id]
//! ```
//! Time slots are the time of reception of the spike, masked by
//! ::num_delay_slots_mask, and neuron IDs are extracted from the spike key by
//! masking with ::incoming_neuron_mask
static uint8_t **spike_counters = NULL;
//! \brief Array of bitfields describing which neurons to deliver spikes to,
//! from which bucket
static bit_field_t *neuron_delay_stage_config = NULL;
//! The number of delay stages. A power of 2.
static uint32_t num_delay_stages = 0;
//! Mask for converting time into the current delay slot
static uint32_t num_delay_slots_mask = 0;
//! Size of each bitfield in ::neuron_delay_stage_config
static uint32_t neuron_bit_field_words = 0;

//! Number of input spikes
static uint32_t n_in_spikes = 0;
//! Number of spikes transferred via queue
static uint32_t n_processed_spikes = 0;
//! Number of spikes sent
static uint32_t n_spikes_sent = 0;
//! Number of spikes added to delay processing
static uint32_t n_spikes_added = 0;

//! \brief An amount of microseconds to back off before starting the timer, in
//! an attempt to avoid overloading the network
static uint32_t timer_offset;

//! \brief The number of clock ticks between processing each neuron at each
//! delay stage
static uint32_t time_between_spikes;

//! The expected current clock tick of timer_1 to wait for
static uint32_t expected_time;

//! Number of times we had to back off because the comms hardware was busy
static uint32_t n_delays = 0;

//! Spin1 API ticks - to know when the timer wraps
extern uint ticks;

//! Used for configuring the timer hardware
static uint32_t timer_period = 0;

//---------------------------------------
// Because we don't want to include string.h or strings.h for memset
//! \brief Sets an array of counters to zero
//!
//! This is basically just bzero()
//!
//! \param[out] counters: The array to zero
//! \param[in] num_items: The size of the array
static inline void zero_spike_counters(
        uint8_t *counters, uint32_t num_items) {
    for (uint32_t i = 0 ; i < num_items ; i++) {
        counters[i] = 0;
    }
}

//! \brief Rounds up to the next power of two
//! \param[in] v: The value to round up
//! \return The minimum power of two that is no smaller than the argument
static inline uint32_t round_to_next_pot(uint32_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

//! \brief Read the configuration region.
//! \param[in] params: The configuration region.
//! \return True if successful
static bool read_parameters(struct delay_parameters *params) {
    log_mini_debug("%u", 9307);  /* log_debug("read_parameters: starting");*/

    key = params->key;
    incoming_key = params->incoming_key;
    incoming_mask = params->incoming_mask;
    incoming_neuron_mask = ~incoming_mask;
    log_mini_debug("%u%08x%08x%08x%08x", 9308, key,  incoming_key,  incoming_mask,  incoming_neuron_mask);
    /* log_debug("\t key = 0x%08x, incoming key = 0x%08x, incoming mask = 0x%08x,incoming key mask = 0x%08x",key, incoming_key, incoming_mask, incoming_neuron_mask);*/


    num_neurons = params->n_atoms;
    neuron_bit_field_words = get_bit_field_size(num_neurons);

    num_delay_stages = params->n_delay_stages;
    timer_offset = params->random_backoff;
    time_between_spikes = params->time_between_spikes * sv->cpu_clk;

    uint32_t num_delay_slots = num_delay_stages * DELAY_STAGE_LENGTH;
    uint32_t num_delay_slots_pot = round_to_next_pot(num_delay_slots);
    num_delay_slots_mask = num_delay_slots_pot - 1;

    log_mini_debug("%u%u%u%u%u%u%08x", 9309, num_neurons,  neuron_bit_field_words, num_delay_stages,  num_delay_slots,  num_delay_slots_pot, num_delay_slots_mask);
    /* log_debug("\t parrot neurons = %u, neuron bit field words = %u, num delay stages = %u, num delay slots = %u (pot = %u), num delay slots mask = %08x",num_neurons, neuron_bit_field_words,num_delay_stages, num_delay_slots, num_delay_slots_pot,num_delay_slots_mask);*/





    log_mini_debug("%u%u%u", 9310, timer_offset,  time_between_spikes);
    /* log_debug("\t random back off = %u, time_between_spikes = %u",timer_offset, time_between_spikes);*/

    // Create array containing a bitfield specifying whether each neuron should
    // emit spikes after each delay stage
    neuron_delay_stage_config =
            spin1_malloc(num_delay_stages * sizeof(bit_field_t));
    if (neuron_delay_stage_config == NULL) {
        log_mini_error("%u%u", 9311, num_delay_stages * sizeof(bit_field_t));
        /* log_error("failed to allocate memory for array of size %u bytes",num_delay_stages * sizeof(bit_field_t));*/
        return false;
    }

    // Loop through delay stages
    for (uint32_t d = 0; d < num_delay_stages; d++) {
        log_mini_debug("%u%u", 9312,  d);  /* log_debug("\t delay stage %u", d);*/

        // Allocate bit-field
        neuron_delay_stage_config[d] =
                spin1_malloc(neuron_bit_field_words * sizeof(uint32_t));
        if (neuron_delay_stage_config[d] == NULL) {
            log_mini_error("%u%u", 9313, neuron_bit_field_words * sizeof(uint32_t));
            /* log_error("failed to allocate memory for bitfield of size %u bytes",neuron_bit_field_words * sizeof(uint32_t));*/
            return false;
        }

        // Copy delay stage configuration bits into delay stage configuration
        // bit-field
        spin1_memcpy(neuron_delay_stage_config[d],
                &params->delay_blocks[d * neuron_bit_field_words],
                neuron_bit_field_words * sizeof(uint32_t));

        for (uint32_t w = 0; w < neuron_bit_field_words; w++) {
            log_mini_debug("%u%u%08x", 9314, w,  neuron_delay_stage_config[d][w]);
            /* log_debug("\t\t delay stage config word %u = %08x",w, neuron_delay_stage_config[d][w]);*/
        }
    }

    // Allocate array of counters for each delay slot
    spike_counters = spin1_malloc(num_delay_slots_pot * sizeof(uint8_t*));
    if (spike_counters == NULL) {
        log_mini_error("%u%u", 9315, num_delay_slots_pot * sizeof(uint8_t*));
        /* log_error("failed to allocate memory for array of size %u bytes",num_delay_slots_pot * sizeof(uint8_t*));*/
        return false;
    }

    for (uint32_t s = 0; s < num_delay_slots_pot; s++) {
        // Allocate an array of counters for each neuron and zero
        spike_counters[s] = spin1_malloc(num_neurons * sizeof(uint8_t));
        if (spike_counters[s] == NULL) {
            log_mini_error("%u%u", 9316, num_neurons * sizeof(uint8_t));
            /* log_error("failed to allocate memory for bitfield of size %u bytes",num_neurons * sizeof(uint8_t));*/
            return false;
        }
        zero_spike_counters(spike_counters[s], num_neurons);
    }

    log_mini_debug("%u", 9317);  /* log_debug("read_parameters: completed successfully");*/
    return true;
}

//! \brief Writes the provenance data
//! \param[out] provenance_region: Where to write the provenance
static void store_provenance_data(address_t provenance_region) {
    log_mini_debug("%u", 9318);  /* log_debug("writing other provenance data");*/
    struct delay_extension_provenance *prov = (void *) provenance_region;

    // store the data into the provenance data region
    prov->n_packets_received = n_in_spikes;
    prov->n_packets_processed = n_processed_spikes;
    prov->n_packets_added = n_spikes_added;
    prov->n_packets_sent = n_spikes_sent;
    prov->n_buffer_overflows = in_spikes_get_n_buffer_overflows();
    prov->n_delays = n_delays;
    log_mini_debug("%u", 9319);  /* log_debug("finished other provenance data");*/
}

//! \brief Read the application configuration
//! \return True if initialisation succeeded.
static bool initialize(void) {
    log_mini_info("%u", 9320);  /* log_info("initialise: started");*/

    // Get the address this core's DTCM data starts at from SRAM
    data_specification_metadata_t *ds_regions =
            data_specification_get_data_address();

    // Read the header
    if (!data_specification_read_header(ds_regions)) {
        return false;
    }

    // Get the timing details and set up the simulation interface
    if (!simulation_initialise(
            data_specification_get_region(SYSTEM, ds_regions),
            APPLICATION_NAME_HASH, &timer_period, &simulation_ticks,
            &infinite_run, &time, SDP, DMA)) {
        return false;
    }
    simulation_set_provenance_function(
            store_provenance_data,
            data_specification_get_region(PROVENANCE_REGION, ds_regions));

    // Get the parameters
    if (!read_parameters(data_specification_get_region(
            DELAY_PARAMS, ds_regions))) {
        return false;
    }

    log_mini_info("%u", 9321);  /* log_info("initialise: completed successfully");*/

    return true;
}

// Callbacks
//! \brief Handles incoming spikes (FIQ)
//!
//! Adds the spikes to the circular buffer handling spikes for later handling by
//! ::spike_process()
//!
//! \param[in] key: the key of the multicast message
//! \param payload: ignored
static void incoming_spike_callback(uint key, uint payload) {
    use(payload);

    log_mini_debug("%u%x", 9322,  key);  /* log_debug("Received spike %x", key);*/
    n_in_spikes++;

    // If there was space to add spike to incoming spike queue
    in_spikes_add_spike(key);
}

//! \brief Gets the neuron ID of the incoming spike
//! \param[in] k: The key
//! \return the neuron ID
static inline index_t key_n(key_t k) {
    return k & incoming_neuron_mask;
}

//! \brief Processes spikes queued by ::incoming_spike_callback()
//!
//! Note that this has to be fairly fast; it is processing with interrupts off.
static inline void spike_process(void) {
    // turn off interrupts as this function is critical for
    // keeping time in sync.
    uint state = spin1_int_disable();

    // Get current time slot of incoming spike counters
    uint32_t current_time_slot = time & num_delay_slots_mask;
    uint8_t *current_time_slot_spike_counters =
            spike_counters[current_time_slot];

    log_mini_debug("%u%u", 9323,  current_time_slot);  /* log_debug("Current time slot %u", current_time_slot);*/

    // While there are any incoming spikes
    spike_t s;
    while (in_spikes_get_next_spike(&s)) {
        n_processed_spikes++;

        if ((s & incoming_mask) == incoming_key) {
            // Mask out neuron ID
            uint32_t neuron_id = key_n(s);
            if (neuron_id < num_neurons) {
                // Increment counter
                current_time_slot_spike_counters[neuron_id]++;
                log_mini_debug("%u%u%u", 9324, neuron_id, current_time_slot_spike_counters[neuron_id]);
                /* log_debug("Incrementing counter %u = %u\n",neuron_id,current_time_slot_spike_counters[neuron_id]);*/

                n_spikes_added++;
            } else {
                log_mini_debug("%u%u", 9325,  neuron_id);  /* log_debug("Invalid neuron ID %u", neuron_id);*/
            }
        } else {
            log_mini_debug("%u%08x", 9326,  s);  /* log_debug("Invalid spike key 0x%08x", s);*/
        }
    }

    // reactivate interrupts as critical section complete
    spin1_mode_restore(state);
}

//! \brief Main timer callback
//! \param[in] timer_count: The current time
//! \param unused1: unused
static void timer_callback(uint timer_count, uint unused1) {
    use(unused1);

    // Process all the spikes from the last timestep
    spike_process();

    time++;

    log_mini_debug("%u%u", 9327,  time);  /* log_debug("Timer tick %u", time);*/

    // If a fixed number of simulation ticks are specified and these have passed
    if (infinite_run != TRUE && time >= simulation_ticks) {
        // handle the pause and resume functionality
        simulation_handle_pause_resume(NULL);

        log_mini_debug("%u%u%u%u%u%u", 9328, time,  n_in_spikes,  n_processed_spikes,  n_spikes_sent, n_spikes_added);
        /* log_debug("Delay extension finished at time %u, %u received spikes, %u processed spikes, %u sent spikes, %u added spikes",time, n_in_spikes, n_processed_spikes, n_spikes_sent,n_spikes_added);*/



        log_mini_debug("%u%u", 9329,  n_delays);  /* log_debug("Delayed %u times", n_delays);*/

        // Subtract 1 from the time so this tick gets done again on the next
        // run
        time--;

        simulation_ready_to_read();
        return;
    }

    // Set the next expected time to wait for between spike sending
    expected_time = sv->cpu_clk * timer_period;

    // Loop through delay stages
    for (uint32_t d = 0; d < num_delay_stages; d++) {
        // If any neurons emit spikes after this delay stage
        bit_field_t delay_stage_config = neuron_delay_stage_config[d];
        if (nonempty_bit_field(delay_stage_config, neuron_bit_field_words)) {
            // Get key mask for this delay stage and it's time slot
            uint32_t delay_stage_delay = (d + 1) * DELAY_STAGE_LENGTH;
            uint32_t delay_stage_time_slot =
                    (time - delay_stage_delay) & num_delay_slots_mask;
            uint8_t *delay_stage_spike_counters =
                    spike_counters[delay_stage_time_slot];

            log_mini_debug("%u%u%u%u", 9330, time,  delay_stage_time_slot,  d);
            /* log_debug("%u: Checking time slot %u for delay stage %u",time, delay_stage_time_slot, d);*/

            // Loop through neurons
            for (uint32_t n = 0; n < num_neurons; n++) {
                // If this neuron emits a spike after this stage
                if (bit_field_test(delay_stage_config, n)) {
                    // Calculate key all spikes coming from this neuron will be
                    // sent with
                    uint32_t spike_key = ((d * num_neurons) + n) + key;

                    if (delay_stage_spike_counters[n] > 0) {
                        log_mini_debug("%u%u%u%u%x", 9331, n,  delay_stage_spike_counters[n],  d, spike_key);
                        /* log_debug("Neuron %u sending %u spikes after delaystage %u with key %x",n, delay_stage_spike_counters[n], d,spike_key);*/


                    }

                    // Loop through counted spikes and send
                    for (uint32_t s = 0; s < delay_stage_spike_counters[n]; s++) {
                        while (!spin1_send_mc_packet(spike_key, 0, NO_PAYLOAD)) {
                            spin1_delay_us(1);
                        }
                        n_spikes_sent++;
                    }
                }

                // Wait until the expected time to send
                while ((ticks == timer_count) && (tc[T1_COUNT] > expected_time)) {
                    // Do Nothing
                    n_delays++;
                }
                expected_time -= time_between_spikes;
            }
        }
    }

    // Zero all counters in current time slot
    uint32_t current_time_slot = time & num_delay_slots_mask;
    zero_spike_counters(spike_counters[current_time_slot], num_neurons);
}

//! Entry point
void c_main(void) {
    if (!initialize()) {
        log_mini_error("%u", 9332);  /* log_error("Error in initialisation - exiting!");*/
        rt_error(RTE_SWERR);
    }

    // Start the time at "-1" so that the first tick will be 0
    time = UINT32_MAX;

    // Initialise the incoming spike buffer
    if (!in_spikes_initialize_spike_buffer(256)) {
        rt_error(RTE_SWERR);
    }

    // Set timer tick (in microseconds)
    log_mini_debug("%u%u", 9333,  timer_period);  /* log_debug("Timer period %u", timer_period);*/
    spin1_set_timer_tick_and_phase(timer_period, timer_offset);

    // Register callbacks
    spin1_callback_on(MC_PACKET_RECEIVED, incoming_spike_callback, MC_PACKET);
    spin1_callback_on(TIMER_TICK, timer_callback, TIMER);

    simulation_run();
}
