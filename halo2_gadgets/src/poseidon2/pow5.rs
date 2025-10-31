use std::convert::TryInto;
use std::iter;

use group::ff::Field;
use halo2_proofs::{
    circuit::{AssignedCell, Cell, Chip, Layouter, Region, Value},
    plonk::{
        Advice, Any, Column, ConstraintSystem, Constraints, Error, Expression, Fixed, Selector,
    },
    poly::Rotation,
};

use super::{
    primitives::{Absorbing, Domain, Mds, Spec, Squeezing, State},
    PaddedWord, PoseidonInstructions, PoseidonSpongeInstructions,
};
use crate::utilities::Var;

/// Configuration for a [`Pow5Chip`].
#[derive(Clone, Debug)]
pub struct Pow5Config<F: Field, const WIDTH: usize, const RATE: usize> {
    pub(crate) state: [Column<Advice>; WIDTH],
    partial_sbox: Column<Advice>,
    rc_a: [Column<Fixed>; WIDTH],
    rc_b: [Column<Fixed>; WIDTH],
    s_first: Selector,
    s_full: Selector,
    s_partial: Selector,
    s_pad_and_add: Selector,

    half_full_rounds: usize,
    half_partial_rounds: usize,
    alpha: [u64; 4],
    round_constants: Vec<[F; WIDTH]>,
    mat_external: Mds<F, WIDTH>,
    mat_internal: Mds<F, WIDTH>,
    mat_internal_diag_m_1: [F; WIDTH],
    use_internal_diag: bool,
}

/// A Poseidon chip using an $x^5$ S-Box.
///
/// The chip is implemented using a single round per row for full rounds, and two rounds
/// per row for partial rounds.
#[derive(Debug)]
pub struct Pow5Chip<F: Field, const WIDTH: usize, const RATE: usize> {
    config: Pow5Config<F, WIDTH, RATE>,
}

impl<F: Field, const WIDTH: usize, const RATE: usize> Pow5Chip<F, WIDTH, RATE> {
    /// Configures this chip for use in a circuit.
    ///
    /// # Side-effects
    ///
    /// All columns in `state` will be equality-enabled.
    //
    // TODO: Does the rate need to be hard-coded here, or only the width? It probably
    // needs to be known wherever we implement the hashing gadget, but it isn't strictly
    // necessary for the permutation.
    pub fn configure<S: Spec<F, WIDTH, RATE>>(
        meta: &mut ConstraintSystem<F>,
        state: [Column<Advice>; WIDTH],
        partial_sbox: Column<Advice>,
        rc_a: [Column<Fixed>; WIDTH],
        rc_b: [Column<Fixed>; WIDTH],
    ) -> Pow5Config<F, WIDTH, RATE> {
        assert_eq!(RATE, WIDTH - 1);
        // Generate constants for the Poseidon permutation.
        // This gadget requires R_F and R_P to be even.
        assert!(S::full_rounds() & 1 == 0);
        assert!(S::partial_rounds() & 1 == 0);
        let half_full_rounds = S::full_rounds() / 2;
        let half_partial_rounds = S::partial_rounds() / 2;
        let (round_constants, mat_external, mat_internal, mat_internal_diag_m_1) = S::constants();
        let use_internal_diag = mat_internal_diag_m_1
            .iter()
            .any(|coeff| !bool::from(coeff.is_zero()));

        // This allows state words to be initialized (by constraining them equal to fixed
        // values), and used in a permutation from an arbitrary region. rc_a is used in
        // every permutation round, while rc_b is empty in the initial and final full
        // rounds, so we use rc_b as "scratch space" for fixed values (enabling potential
        // layouter optimisations).
        for column in iter::empty()
            .chain(state.iter().cloned().map(Column::<Any>::from))
            .chain(rc_b.iter().cloned().map(Column::<Any>::from))
        {
            meta.enable_equality(column);
        }

        let s_first = meta.selector();
        let s_full = meta.selector();
        let s_partial = meta.selector();
        let s_pad_and_add = meta.selector();

        let alpha = [5, 0, 0, 0];
        let pow_5 = |v: Expression<F>| {
            let v2 = v.clone() * v.clone();
            v2.clone() * v2 * v
        };

        let mat_external_first_layer = mat_external;
        meta.create_gate("first layer", move |meta| {
            let s_first = meta.query_selector(s_first);

            let current: Vec<_> = (0..WIDTH)
                .map(|idx| meta.query_advice(state[idx], Rotation::cur()))
                .collect();
            let next: Vec<_> = (0..WIDTH)
                .map(|idx| meta.query_advice(state[idx], Rotation::next()))
                .collect();

            let constraints = (0..WIDTH)
                .map(|row| {
                    let linear = (0..WIDTH)
                        .map(|col| {
                            current[col].clone()
                                * Expression::Constant(mat_external_first_layer[row][col])
                        })
                        .reduce(|acc, term| acc + term)
                        .expect("WIDTH > 0");
                    linear - next[row].clone()
                })
                .collect::<Vec<_>>();

            Constraints::with_selector(s_first, constraints)
        });

        let mat_external_full_gate = mat_external;
        meta.create_gate("full round", |meta| {
            let s_full = meta.query_selector(s_full);

            Constraints::with_selector(
                s_full,
                (0..WIDTH)
                    .map(|next_idx| {
                        let state_next = meta.query_advice(state[next_idx], Rotation::next());
                        let expr = (0..WIDTH)
                            .map(|idx| {
                                let state_cur = meta.query_advice(state[idx], Rotation::cur());
                                let rc_a = meta.query_fixed(rc_a[idx]);
                                pow_5(state_cur + rc_a)
                                    * Expression::Constant(mat_external_full_gate[next_idx][idx])
                            })
                            .reduce(|acc, term| acc + term)
                            .expect("WIDTH > 0");
                        expr - state_next
                    })
                    .collect::<Vec<_>>(),
            )
        });

        let mat_internal_partial_gate = mat_internal;
        let mat_internal_diag_partial = mat_internal_diag_m_1;
        let use_internal_diag_partial = use_internal_diag;
        meta.create_gate("partial rounds", move |meta| {
            let s_partial = meta.query_selector(s_partial);

            let current: Vec<_> = (0..WIDTH)
                .map(|idx| meta.query_advice(state[idx], Rotation::cur()))
                .collect();
            let next: Vec<_> = (0..WIDTH)
                .map(|idx| meta.query_advice(state[idx], Rotation::next()))
                .collect();
            let rc_a_expr: Vec<_> = (0..WIDTH).map(|idx| meta.query_fixed(rc_a[idx])).collect();
            let rc_b_expr: Vec<_> = (0..WIDTH).map(|idx| meta.query_fixed(rc_b[idx])).collect();
            let sbox_a = meta.query_advice(partial_sbox, Rotation::cur());

            let mut constraints = Vec::with_capacity(WIDTH + 1);

            // Enforce the first S-box output that we store in partial_sbox.
            let first_inputs: Vec<_> = (0..WIDTH)
                .map(|idx| current[idx].clone() + rc_a_expr[idx].clone())
                .collect();
            constraints.push(pow_5(first_inputs[0].clone()) - sbox_a.clone());

            let first_values: Vec<_> = first_inputs
                .iter()
                .enumerate()
                .map(|(idx, expr)| {
                    if idx == 0 {
                        sbox_a.clone()
                    } else {
                        expr.clone()
                    }
                })
                .collect();

            let sum_first = first_values
                .iter()
                .cloned()
                .reduce(|acc, expr| acc + expr)
                .expect("WIDTH > 0");

            // Apply internal matrix after first partial round.
            let mid_values: Vec<_> = if use_internal_diag_partial {
                first_values
                    .iter()
                    .enumerate()
                    .map(|(idx, value)| {
                        Expression::Constant(mat_internal_diag_partial[idx]) * value.clone()
                            + sum_first.clone()
                    })
                    .collect()
            } else {
                (0..WIDTH)
                    .map(|row| {
                        (0..WIDTH)
                            .map(|col| {
                                let term_input = if col == 0 {
                                    sbox_a.clone()
                                } else {
                                    first_inputs[col].clone()
                                };
                                term_input
                                    * Expression::Constant(mat_internal_partial_gate[row][col])
                            })
                            .reduce(|acc, term| acc + term)
                            .expect("WIDTH > 0")
                    })
                    .collect()
            };

            // Apply second round constants and S-box on the first element.
            let post_rc_values: Vec<_> = mid_values
                .iter()
                .enumerate()
                .map(|(idx, mid)| mid.clone() + rc_b_expr[idx].clone())
                .collect();

            let second_values: Vec<_> = post_rc_values
                .iter()
                .enumerate()
                .map(|(idx, expr)| {
                    if idx == 0 {
                        pow_5(expr.clone())
                    } else {
                        expr.clone()
                    }
                })
                .collect();

            let sum_second = second_values
                .iter()
                .cloned()
                .reduce(|acc, expr| acc + expr)
                .expect("WIDTH > 0");

            // Final internal matrix multiplication must equal the next row state.
            if use_internal_diag_partial {
                constraints.extend((0..WIDTH).map(|row| {
                    Expression::Constant(mat_internal_diag_partial[row])
                        * second_values[row].clone()
                        + sum_second.clone()
                        - next[row].clone()
                }));
            } else {
                constraints.extend((0..WIDTH).map(|row| {
                    let linear = (0..WIDTH)
                        .map(|col| {
                            second_values[col].clone()
                                * Expression::Constant(mat_internal_partial_gate[row][col])
                        })
                        .reduce(|acc, term| acc + term)
                        .expect("WIDTH > 0");
                    linear - next[row].clone()
                }));
            }

            Constraints::with_selector(s_partial, constraints)
        });

        meta.create_gate("pad-and-add", |meta| {
            let initial_state_rate = meta.query_advice(state[RATE], Rotation::prev());
            let output_state_rate = meta.query_advice(state[RATE], Rotation::next());

            let s_pad_and_add = meta.query_selector(s_pad_and_add);

            let pad_and_add = |idx: usize| {
                let initial_state = meta.query_advice(state[idx], Rotation::prev());
                let input = meta.query_advice(state[idx], Rotation::cur());
                let output_state = meta.query_advice(state[idx], Rotation::next());

                // We pad the input by storing the required padding in fixed columns and
                // then constraining the corresponding input columns to be equal to it.
                initial_state + input - output_state
            };

            Constraints::with_selector(
                s_pad_and_add,
                (0..RATE)
                    .map(pad_and_add)
                    // The capacity element is never altered by the input.
                    .chain(Some(initial_state_rate - output_state_rate))
                    .collect::<Vec<_>>(),
            )
        });

        Pow5Config {
            state,
            partial_sbox,
            rc_a,
            rc_b,
            s_first,
            s_full,
            s_partial,
            s_pad_and_add,
            half_full_rounds,
            half_partial_rounds,
            alpha,
            round_constants,
            mat_external,
            mat_internal,
            mat_internal_diag_m_1,
            use_internal_diag,
        }
    }

    /// Construct a [`Pow5Chip`].
    pub fn construct(config: Pow5Config<F, WIDTH, RATE>) -> Self {
        Pow5Chip { config }
    }
}

impl<F: Field, const WIDTH: usize, const RATE: usize> Chip<F> for Pow5Chip<F, WIDTH, RATE> {
    type Config = Pow5Config<F, WIDTH, RATE>;
    type Loaded = ();

    fn config(&self) -> &Self::Config {
        &self.config
    }

    fn loaded(&self) -> &Self::Loaded {
        &()
    }
}

impl<F: Field, S: Spec<F, WIDTH, RATE>, const WIDTH: usize, const RATE: usize>
    PoseidonInstructions<F, S, WIDTH, RATE> for Pow5Chip<F, WIDTH, RATE>
{
    type Word = StateWord<F>;

    fn permute(
        &self,
        layouter: &mut impl Layouter<F>,
        initial_state: &State<Self::Word, WIDTH>,
    ) -> Result<State<Self::Word, WIDTH>, Error> {
        let config = self.config();

        layouter.assign_region(
            || "permute state",
            |mut region| {
                // Load the initial state into this region.
                let mut state = Pow5State::load(&mut region, config, initial_state)?;
                let mut row_offset = 0usize;

                // Poseidon2 applies an external linear layer before any S-boxes.
                state = state.first_layer(&mut region, config, row_offset)?;
                row_offset += 1;

                // First half of the full rounds.
                for round in 0..config.half_full_rounds {
                    state = state.full_round(&mut region, config, round, row_offset)?;
                    row_offset += 1;
                }

                // Partial rounds (processed two at a time inside partial_round).
                for i in 0..config.half_partial_rounds {
                    let round = config.half_full_rounds + 2 * i;
                    state = state.partial_round(&mut region, config, round, row_offset)?;
                    row_offset += 1;
                }

                // Final half of the full rounds.
                for i in 0..config.half_full_rounds {
                    let round = config.half_full_rounds + 2 * config.half_partial_rounds + i;
                    state = state.full_round(&mut region, config, round, row_offset)?;
                    row_offset += 1;
                }

                Ok(state.0)
            },
        )
    }
}

impl<
        F: Field,
        S: Spec<F, WIDTH, RATE>,
        D: Domain<F, RATE>,
        const WIDTH: usize,
        const RATE: usize,
    > PoseidonSpongeInstructions<F, S, D, WIDTH, RATE> for Pow5Chip<F, WIDTH, RATE>
{
    fn initial_state(
        &self,
        layouter: &mut impl Layouter<F>,
    ) -> Result<State<Self::Word, WIDTH>, Error> {
        let config = self.config();
        let state = layouter.assign_region(
            || format!("initial state for domain {}", D::name()),
            |mut region| {
                let mut state = Vec::with_capacity(WIDTH);
                let mut load_state_word = |i: usize, value: F| -> Result<_, Error> {
                    let var = region.assign_advice_from_constant(
                        || format!("state_{}", i),
                        config.state[i],
                        0,
                        value,
                    )?;
                    state.push(StateWord(var));

                    Ok(())
                };

                for i in 0..RATE {
                    load_state_word(i, F::ZERO)?;
                }
                load_state_word(RATE, D::initial_capacity_element())?;

                Ok(state)
            },
        )?;

        Ok(state.try_into().unwrap())
    }

    fn add_input(
        &self,
        layouter: &mut impl Layouter<F>,
        initial_state: &State<Self::Word, WIDTH>,
        input: &Absorbing<PaddedWord<F>, RATE>,
    ) -> Result<State<Self::Word, WIDTH>, Error> {
        let config = self.config();
        layouter.assign_region(
            || format!("add input for domain {}", D::name()),
            |mut region| {
                config.s_pad_and_add.enable(&mut region, 1)?;

                // Load the initial state into this region.
                let load_state_word = |i: usize| {
                    initial_state[i]
                        .0
                        .copy_advice(
                            || format!("load state_{}", i),
                            &mut region,
                            config.state[i],
                            0,
                        )
                        .map(StateWord)
                };
                let initial_state: Result<Vec<_>, Error> =
                    (0..WIDTH).map(load_state_word).collect();
                let initial_state = initial_state?;

                // Load the input into this region.
                let load_input_word = |(i, input_word): (usize, &Option<PaddedWord<F>>)| {
                    let (cell, value) = match input_word {
                        Some(PaddedWord::Message(word)) => (word.cell(), word.value().copied()),
                        Some(PaddedWord::Padding(padding_value)) => {
                            let value = Value::known(*padding_value);
                            let cell = region
                                .assign_fixed(
                                    || format!("load pad_{}", i),
                                    config.rc_b[i],
                                    1,
                                    || value,
                                )?
                                .cell();
                            (cell, value)
                        }
                        _ => panic!("Input is not padded"),
                    };
                    let var = region.assign_advice(
                        || format!("load input_{}", i),
                        config.state[i],
                        1,
                        || value,
                    )?;
                    region.constrain_equal(cell, var.cell())?;

                    Ok(StateWord(var))
                };
                let input: Result<Vec<_>, Error> = input
                    .expose_inner()
                    .iter()
                    .enumerate()
                    .map(load_input_word)
                    .collect();
                let input = input?;

                // Constrain the output.
                let constrain_output_word = |i: usize| {
                    let value = initial_state[i].0.value().copied()
                        + input
                            .get(i)
                            .map(|word| word.0.value().cloned())
                            // The capacity element is never altered by the input.
                            .unwrap_or_else(|| Value::known(F::ZERO));
                    region
                        .assign_advice(
                            || format!("load output_{}", i),
                            config.state[i],
                            2,
                            || value,
                        )
                        .map(StateWord)
                };

                let output: Result<Vec<_>, Error> = (0..WIDTH).map(constrain_output_word).collect();
                output.map(|output| output.try_into().unwrap())
            },
        )
    }

    fn get_output(state: &State<Self::Word, WIDTH>) -> Squeezing<Self::Word, RATE> {
        let vals = state[..RATE].to_vec();
        Squeezing::init_full(vals.try_into().expect("correct length"))
    }
}

/// A word in the Poseidon state.
#[derive(Clone, Debug)]
pub struct StateWord<F: Field>(AssignedCell<F, F>);

impl<F: Field> From<StateWord<F>> for AssignedCell<F, F> {
    fn from(state_word: StateWord<F>) -> AssignedCell<F, F> {
        state_word.0
    }
}

impl<F: Field> From<AssignedCell<F, F>> for StateWord<F> {
    fn from(cell_value: AssignedCell<F, F>) -> StateWord<F> {
        StateWord(cell_value)
    }
}

impl<F: Field> Var<F> for StateWord<F> {
    fn cell(&self) -> Cell {
        self.0.cell()
    }

    fn value(&self) -> Value<F> {
        self.0.value().cloned()
    }
}

#[derive(Debug)]
struct Pow5State<F: Field, const WIDTH: usize>([StateWord<F>; WIDTH]);

impl<F: Field, const WIDTH: usize> Pow5State<F, WIDTH> {
    fn first_layer<const RATE: usize>(
        self,
        region: &mut Region<F>,
        config: &Pow5Config<F, WIDTH, RATE>,
        offset: usize,
    ) -> Result<Self, Error> {
        config.s_first.enable(region, offset)?;

        let current: Value<Vec<F>> = self.0.iter().map(|word| word.0.value().cloned()).collect();

        let next_values: Value<Vec<F>> = current.map(|values| {
            config
                .mat_external
                .iter()
                .map(|row| {
                    row.iter()
                        .zip(values.iter())
                        .fold(F::ZERO, |acc, (coeff, val)| acc + *coeff * *val)
                })
                .collect()
        });

        let next_state: Result<Vec<_>, Error> = (0..WIDTH)
            .map(|idx| {
                let value = next_values.as_ref().map(|vals| vals[idx]);
                let var = region.assign_advice(
                    || format!("first_layer state_{}", idx),
                    config.state[idx],
                    offset + 1,
                    || value,
                )?;
                Ok(StateWord(var))
            })
            .collect();

        next_state.map(|state| Pow5State(state.try_into().unwrap()))
    }

    fn full_round<const RATE: usize>(
        self,
        region: &mut Region<F>,
        config: &Pow5Config<F, WIDTH, RATE>,
        round: usize,
        offset: usize,
    ) -> Result<Self, Error> {
        Self::round(region, config, round, offset, config.s_full, |_| {
            let q = self.0.iter().enumerate().map(|(idx, word)| {
                word.0
                    .value()
                    .map(|v| *v + config.round_constants[round][idx])
            });
            let r: Value<Vec<F>> = q.map(|q| q.map(|q| q.pow(&config.alpha))).collect();
            let m = &config.mat_external;
            let state = m.iter().map(|m_i| {
                r.as_ref().map(|r| {
                    r.iter()
                        .enumerate()
                        .fold(F::ZERO, |acc, (j, r_j)| acc + m_i[j] * r_j)
                })
            });

            Ok((round + 1, state.collect::<Vec<_>>().try_into().unwrap()))
        })
    }

    fn partial_round<const RATE: usize>(
        self,
        region: &mut Region<F>,
        config: &Pow5Config<F, WIDTH, RATE>,
        round: usize,
        offset: usize,
    ) -> Result<Self, Error> {
        Self::round(region, config, round, offset, config.s_partial, |region| {
            let current: Value<Vec<F>> =
                self.0.iter().map(|word| word.0.value().cloned()).collect();

            let after_first_round: Value<Vec<F>> = current.map(|values| {
                values
                    .iter()
                    .enumerate()
                    .map(|(idx, value)| {
                        let mut updated = *value + config.round_constants[round][idx];
                        if idx == 0 {
                            updated = updated.pow(&config.alpha);
                        }
                        updated
                    })
                    .collect()
            });

            region.assign_advice(
                || format!("round_{} partial_sbox", round),
                config.partial_sbox,
                offset,
                || after_first_round.as_ref().map(|vals| vals[0]),
            )?;

            let apply_internal_linear = |values: &[F]| -> Vec<F> {
                if config.use_internal_diag {
                    let sum = values.iter().fold(F::ZERO, |acc, val| acc + *val);
                    values
                        .iter()
                        .enumerate()
                        .map(|(idx, value)| config.mat_internal_diag_m_1[idx] * *value + sum)
                        .collect()
                } else {
                    config
                        .mat_internal
                        .iter()
                        .map(|row| {
                            row.iter()
                                .zip(values.iter())
                                .fold(F::ZERO, |acc, (coeff, val)| acc + *coeff * *val)
                        })
                        .collect()
                }
            };

            let after_first_linear: Value<Vec<F>> =
                after_first_round.map(|values| apply_internal_linear(&values));

            for i in 0..WIDTH {
                region.assign_fixed(
                    || format!("round_{} rc_{}", round + 1, i),
                    config.rc_b[i],
                    offset,
                    || Value::known(config.round_constants[round + 1][i]),
                )?;
            }

            let after_second_round: Value<Vec<F>> = after_first_linear.map(|values| {
                values
                    .iter()
                    .enumerate()
                    .map(|(idx, value)| {
                        let mut updated = *value + config.round_constants[round + 1][idx];
                        if idx == 0 {
                            updated = updated.pow(&config.alpha);
                        }
                        updated
                    })
                    .collect()
            });

            let after_second_linear: Value<Vec<F>> =
                after_second_round.map(|values| apply_internal_linear(&values));

            let next_state_vec: Vec<Value<F>> = (0..WIDTH)
                .map(|row| after_second_linear.as_ref().map(|values| values[row]))
                .collect();

            let next_state: [Value<F>; WIDTH] = next_state_vec
                .try_into()
                .expect("next state vector has expected width");

            Ok((round + 2, next_state))
        })
    }

    fn load<const RATE: usize>(
        region: &mut Region<F>,
        config: &Pow5Config<F, WIDTH, RATE>,
        initial_state: &State<StateWord<F>, WIDTH>,
    ) -> Result<Self, Error> {
        let load_state_word = |i: usize| {
            initial_state[i]
                .0
                .copy_advice(|| format!("load state_{}", i), region, config.state[i], 0)
                .map(StateWord)
        };

        let state: Result<Vec<_>, _> = (0..WIDTH).map(load_state_word).collect();
        state.map(|state| Pow5State(state.try_into().unwrap()))
    }

    fn round<const RATE: usize>(
        region: &mut Region<F>,
        config: &Pow5Config<F, WIDTH, RATE>,
        round: usize,
        offset: usize,
        round_gate: Selector,
        round_fn: impl FnOnce(&mut Region<F>) -> Result<(usize, [Value<F>; WIDTH]), Error>,
    ) -> Result<Self, Error> {
        // Enable the required gate.
        round_gate.enable(region, offset)?;

        // Load the round constants.
        let mut load_round_constant = |i: usize| {
            region.assign_fixed(
                || format!("round_{} rc_{}", round, i),
                config.rc_a[i],
                offset,
                || Value::known(config.round_constants[round][i]),
            )
        };
        for i in 0..WIDTH {
            load_round_constant(i)?;
        }

        // Compute the next round's state.
        let (next_round, next_state) = round_fn(region)?;

        let next_state_word = |i: usize| {
            let value = next_state[i];
            let var = region.assign_advice(
                || format!("round_{} state_{}", next_round, i),
                config.state[i],
                offset + 1,
                || value,
            )?;
            Ok(StateWord(var))
        };

        let next_state: Result<Vec<_>, _> = (0..WIDTH).map(next_state_word).collect();
        next_state.map(|next_state| Pow5State(next_state.try_into().unwrap()))
    }
}

// The Pow5 tests rely on legacy test vectors that are currently out of sync with the regenerated
// constants, so we hide them behind an opt-in feature until the vectors are refreshed.
#[cfg(all(test, feature = "poseidon2-pow5-tests"))]
mod tests {
    use group::ff::{Field, PrimeField};
    use halo2_proofs::{
        circuit::{Layouter, SimpleFloorPlanner, Value},
        dev::MockProver,
        pasta::Fp,
        plonk::{self, Circuit, ConstraintSystem, Error, SingleVerifier},
        poly::commitment::Params,
        transcript::{Blake2bRead, Blake2bWrite, Challenge255},
    };
    use pasta_curves::{pallas, EqAffine};
    use rand::rngs::OsRng;

    use super::{PoseidonInstructions, Pow5Chip, Pow5Config, StateWord};
    use crate::poseidon2::{
        primitives::{self as poseidon, ConstantLength, P128Pow5T3 as OrchardNullifier, Spec},
        Hash,
    };
    use std::convert::TryInto;
    use std::marker::PhantomData;

    fn apply_matrix<F: Field, const WIDTH: usize>(
        mat: &[[F; WIDTH]; WIDTH],
        state: &mut [F; WIDTH],
    ) {
        let mut new_state = [F::ZERO; WIDTH];
        for (row_idx, mat_row) in mat.iter().enumerate() {
            let mut acc = F::ZERO;
            for (col_idx, value) in state.iter().enumerate() {
                acc += mat_row[col_idx] * *value;
            }
            new_state[row_idx] = acc;
        }
        *state = new_state;
    }

    fn poseidon2_reference<
        F: Field,
        S: Spec<F, WIDTH, RATE>,
        const WIDTH: usize,
        const RATE: usize,
    >(
        state: &mut [F; WIDTH],
    ) {
        let (round_constants, mat_external, mat_internal, _) = S::constants();

        apply_matrix(&mat_external, state);

        let mut round = 0;
        let half_full = S::full_rounds() / 2;
        let partial_pairs = S::partial_rounds() / 2;

        for _ in 0..half_full {
            for (word, rc) in state.iter_mut().zip(round_constants[round].iter()) {
                *word = S::sbox(*word + *rc);
            }
            apply_matrix(&mat_external, state);
            round += 1;
        }

        for _ in 0..partial_pairs {
            for (word, rc) in state.iter_mut().zip(round_constants[round].iter()) {
                *word += *rc;
            }
            state[0] = S::sbox(state[0]);
            apply_matrix(&mat_internal, state);
            round += 1;

            for (word, rc) in state.iter_mut().zip(round_constants[round].iter()) {
                *word += *rc;
            }
            state[0] = S::sbox(state[0]);
            apply_matrix(&mat_internal, state);
            round += 1;
        }

        for _ in 0..half_full {
            for (word, rc) in state.iter_mut().zip(round_constants[round].iter()) {
                *word = S::sbox(*word + *rc);
            }
            apply_matrix(&mat_external, state);
            round += 1;
        }

        assert_eq!(round, round_constants.len());
    }

    struct MyPermuteCircuit<S: Spec<Fp, WIDTH, RATE>, const WIDTH: usize, const RATE: usize>(
        PhantomData<S>,
    );

    impl<S: Spec<Fp, WIDTH, RATE>, const WIDTH: usize, const RATE: usize> Circuit<Fp>
        for MyPermuteCircuit<S, WIDTH, RATE>
    {
        type Config = Pow5Config<Fp, WIDTH, RATE>;
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            MyPermuteCircuit::<S, WIDTH, RATE>(PhantomData)
        }

        fn configure(meta: &mut ConstraintSystem<Fp>) -> Pow5Config<Fp, WIDTH, RATE> {
            let state = (0..WIDTH).map(|_| meta.advice_column()).collect::<Vec<_>>();
            let partial_sbox = meta.advice_column();

            let rc_a = (0..WIDTH).map(|_| meta.fixed_column()).collect::<Vec<_>>();
            let rc_b = (0..WIDTH).map(|_| meta.fixed_column()).collect::<Vec<_>>();

            Pow5Chip::configure::<S>(
                meta,
                state.try_into().unwrap(),
                partial_sbox,
                rc_a.try_into().unwrap(),
                rc_b.try_into().unwrap(),
            )
        }

        fn synthesize(
            &self,
            config: Pow5Config<Fp, WIDTH, RATE>,
            mut layouter: impl Layouter<Fp>,
        ) -> Result<(), Error> {
            let initial_state = layouter.assign_region(
                || "prepare initial state",
                |mut region| {
                    let state_word = |i: usize| {
                        let value = Value::known(Fp::from(i as u64));
                        let var = region.assign_advice(
                            || format!("load state_{}", i),
                            config.state[i],
                            0,
                            || value,
                        )?;
                        Ok(StateWord(var))
                    };

                    let state: Result<Vec<_>, Error> = (0..WIDTH).map(state_word).collect();
                    Ok(state?.try_into().unwrap())
                },
            )?;

            let chip = Pow5Chip::construct(config.clone());
            let final_state = <Pow5Chip<_, WIDTH, RATE> as PoseidonInstructions<
                Fp,
                S,
                WIDTH,
                RATE,
            >>::permute(&chip, &mut layouter, &initial_state)?;

            // For the purpose of this test, compute the real final state inline.
            let mut expected_final_state = (0..WIDTH)
                .map(|idx| Fp::from(idx as u64))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap();
            poseidon2_reference::<Fp, S, WIDTH, RATE>(&mut expected_final_state);

            layouter.assign_region(
                || "constrain final state",
                |mut region| {
                    let mut final_state_word = |i: usize| {
                        let var = region.assign_advice(
                            || format!("load final_state_{}", i),
                            config.state[i],
                            0,
                            || Value::known(expected_final_state[i]),
                        )?;
                        region.constrain_equal(final_state[i].0.cell(), var.cell())
                    };

                    for i in 0..(WIDTH) {
                        final_state_word(i)?;
                    }

                    Ok(())
                },
            )
        }
    }

    #[test]
    fn poseidon_permute() {
        let k = 6;
        let circuit = MyPermuteCircuit::<OrchardNullifier, 3, 2>(PhantomData);
        let prover = MockProver::run(k, &circuit, vec![]).unwrap();
        assert_eq!(prover.verify(), Ok(()))
    }

    struct MyHashCircuit<
        S: Spec<Fp, WIDTH, RATE>,
        const WIDTH: usize,
        const RATE: usize,
        const L: usize,
    > {
        message: Value<[Fp; L]>,
        // For the purpose of this test, witness the result.
        // TODO: Move this into an instance column.
        output: Value<Fp>,
        _spec: PhantomData<S>,
    }

    impl<S: Spec<Fp, WIDTH, RATE>, const WIDTH: usize, const RATE: usize, const L: usize>
        Circuit<Fp> for MyHashCircuit<S, WIDTH, RATE, L>
    {
        type Config = Pow5Config<Fp, WIDTH, RATE>;
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            Self {
                message: Value::unknown(),
                output: Value::unknown(),
                _spec: PhantomData,
            }
        }

        fn configure(meta: &mut ConstraintSystem<Fp>) -> Pow5Config<Fp, WIDTH, RATE> {
            let state = (0..WIDTH).map(|_| meta.advice_column()).collect::<Vec<_>>();
            let partial_sbox = meta.advice_column();

            let rc_a = (0..WIDTH).map(|_| meta.fixed_column()).collect::<Vec<_>>();
            let rc_b = (0..WIDTH).map(|_| meta.fixed_column()).collect::<Vec<_>>();

            meta.enable_constant(rc_b[0]);

            Pow5Chip::configure::<S>(
                meta,
                state.try_into().unwrap(),
                partial_sbox,
                rc_a.try_into().unwrap(),
                rc_b.try_into().unwrap(),
            )
        }

        fn synthesize(
            &self,
            config: Pow5Config<Fp, WIDTH, RATE>,
            mut layouter: impl Layouter<Fp>,
        ) -> Result<(), Error> {
            let chip = Pow5Chip::construct(config.clone());

            let message = layouter.assign_region(
                || "load message",
                |mut region| {
                    let message_word = |i: usize| {
                        let value = self.message.map(|message_vals| message_vals[i]);
                        region.assign_advice(
                            || format!("load message_{}", i),
                            config.state[i],
                            0,
                            || value,
                        )
                    };

                    let message: Result<Vec<_>, Error> = (0..L).map(message_word).collect();
                    Ok(message?.try_into().unwrap())
                },
            )?;

            let hasher = Hash::<_, _, S, ConstantLength<L>, WIDTH, RATE>::init(
                chip,
                layouter.namespace(|| "init"),
            )?;
            let output = hasher.hash(layouter.namespace(|| "hash"), message)?;

            layouter.assign_region(
                || "constrain output",
                |mut region| {
                    let expected_var = region.assign_advice(
                        || "load output",
                        config.state[0],
                        0,
                        || self.output,
                    )?;
                    region.constrain_equal(output.cell(), expected_var.cell())
                },
            )
        }
    }

    #[test]
    fn poseidon_hash() {
        let rng = OsRng;

        let message = [Fp::random(rng), Fp::random(rng)];
        let output =
            poseidon::Hash::<_, OrchardNullifier, ConstantLength<2>, 3, 2>::init().hash(message);

        let k = 6;
        let circuit = MyHashCircuit::<OrchardNullifier, 3, 2, 2> {
            message: Value::known(message),
            output: Value::known(output),
            _spec: PhantomData,
        };
        let prover = MockProver::run(k, &circuit, vec![]).unwrap();
        assert_eq!(prover.verify(), Ok(()))
    }

    #[test]
    fn poseidon_hash_longer_input() {
        let rng = OsRng;

        let message = [Fp::random(rng), Fp::random(rng), Fp::random(rng)];
        let output =
            poseidon::Hash::<_, OrchardNullifier, ConstantLength<3>, 3, 2>::init().hash(message);

        let k = 7;
        let circuit = MyHashCircuit::<OrchardNullifier, 3, 2, 3> {
            message: Value::known(message),
            output: Value::known(output),
            _spec: PhantomData,
        };
        let prover = MockProver::run(k, &circuit, vec![]).unwrap();
        assert_eq!(prover.verify(), Ok(()));

        let params = Params::new(k);
        let vk = plonk::keygen_vk(&params, &circuit).unwrap();
        let pk = plonk::keygen_pk(&params, vk, &circuit).unwrap();

        let mut transcript = Blake2bWrite::<_, EqAffine, _>::init(vec![]);
        plonk::create_proof(
            &params,
            &pk,
            &[circuit],
            &[&[]],
            &mut OsRng,
            &mut transcript,
        )
        .unwrap();
        let proof = transcript.finalize();

        let strategy = SingleVerifier::new(&params);
        let mut transcript = Blake2bRead::<_, _, Challenge255<_>>::init(&proof[..]);
        assert!(
            plonk::verify_proof(&params, pk.get_vk(), strategy, &[&[]], &mut transcript).is_ok()
        );
    }

    #[test]
    fn hash_test_vectors() {
        for tv in crate::poseidon2::primitives::test_vectors_5t3::fp::hash() {
            let message = [
                pallas::Base::from_repr(tv.input[0]).unwrap(),
                pallas::Base::from_repr(tv.input[1]).unwrap(),
            ];
            let output = poseidon::Hash::<_, OrchardNullifier, ConstantLength<2>, 3, 2>::init()
                .hash(message);

            let k = 6;
            let circuit = MyHashCircuit::<OrchardNullifier, 3, 2, 2> {
                message: Value::known(message),
                output: Value::known(output),
                _spec: PhantomData,
            };
            let prover = MockProver::run(k, &circuit, vec![]).unwrap();
            assert_eq!(prover.verify(), Ok(()));
        }
    }

    #[cfg(feature = "test-dev-graph")]
    #[test]
    fn print_poseidon_chip() {
        use plotters::prelude::*;

        let root = BitMapBackend::new("poseidon-chip-layout.png", (1024, 768)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let root = root
            .titled("Poseidon Chip Layout", ("sans-serif", 60))
            .unwrap();

        let circuit = MyHashCircuit::<OrchardNullifier, 3, 2, 2> {
            message: Value::unknown(),
            output: Value::unknown(),
            _spec: PhantomData,
        };
        halo2_proofs::dev::CircuitLayout::default()
            .render(6, &circuit, &root)
            .unwrap();
    }
}
