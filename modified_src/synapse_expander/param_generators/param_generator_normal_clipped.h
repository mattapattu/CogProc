// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/synapse_expander/param_generators/param_generator_normal_clipped.h

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

/**
 * \file
 * \brief Normally distributed random, redrawn if out of boundary, parameter
 *        generator implementation
 */
#include <stdfix.h>
#include <spin1_api.h>
#include <stdfix-full-iso.h>
#include <normal.h>
#include <synapse_expander/rng.h>
#include <synapse_expander/generator_types.h>

/**
 * \brief The parameters that can be copied in from SDRAM
 */
struct normal_clipped_params {
    accum mu;
    accum sigma;
    accum low;
    accum high;
};

/**
 * \brief The data structure to be passed around for this generator.  This
 *        includes the parameters and an RNG.
 */
struct param_generator_normal_clipped {
    struct normal_clipped_params params;
    rng_t rng;
};

/**
 * \brief How to initialise the clipped normal RNG parameter generator
 * \param[in,out] region: Region to read setup from.  Should be updated
 *                        to position just after parameters after calling.
 * \return A data item to be passed in to other functions later on
 */
static void *param_generator_normal_clipped_initialize(address_t *region) {
    // Allocate memory for the data
    struct param_generator_normal_clipped *obj =
            spin1_malloc(sizeof(struct param_generator_normal_clipped));
    struct normal_clipped_params *params_sdram = (void *) *region;

    // Copy the parameters in
    obj->params = *params_sdram++;
    *region = (void *) params_sdram;

    log_mini_debug("%u%k%k%k%k", 9362, obj->params.mu,  obj->params.sigma,  obj->params.low, obj->params.high);
    /* log_debug("normal clipped mu = %k, sigma = %k, low = %k, high = %k",obj->params.mu, obj->params.sigma, obj->params.low,obj->params.high);*/


    // Initialise the RNG for this generator
    obj->rng = rng_init(region);
    return obj;
}

/**
 * \brief How to free any data for the clipped normal RNG parameter generator
 * \param[in] generator: The generator to free
 */
static void param_generator_normal_clipped_free(void *generator) {
    struct param_generator_normal_clipped *obj = generator;
    rng_free(obj->rng);
    sark_free(generator);
}

/**
 * \brief How to generate values with the clipped normal RNG parameter generator
 * \param[in] generator: The generator to use to generate values
 * \param[in] n_indices: The number of values to generate
 * \param[in] pre_neuron_index: The index of the neuron in the pre-population
 *                              being generated
 * \param[in] indices: The \p n_indices post-neuron indices for each connection
 * \param[out] values: An array into which to place the values; will be
 *                     \p n_indices in size
 */
static void param_generator_normal_clipped_generate(
        void *generator, uint32_t n_indices, uint32_t pre_neuron_index,
        uint16_t *indices, accum *values) {
    use(pre_neuron_index);
    use(indices);

    // For each index, generate a normally distributed random value, redrawing
    // if outside the given range
    struct param_generator_normal_clipped *obj = generator;
    for (uint32_t i = 0; i < n_indices; i++) {
        do {
            accum value = rng_normal(obj->rng);
            values[i] = obj->params.mu + (value * obj->params.sigma);
        } while (values[i] < obj->params.low || values[i] > obj->params.high);
    }
}
