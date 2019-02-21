
#include "coulomb.hpp"
#include "communication.hpp"                               // bcast functions
#include "electrostatics_magnetostatics/maggs.hpp"         // maggs_init
// #include "electrostatics_magnetostatics/mdlc_correction.hpp"
#include "electrostatics_magnetostatics/p3m_gpu.hpp"  // p3m_gpu_init
#include "errorhandling.hpp"                          // runtime_error
#include "initialize.hpp"                             // on_ghost_flags_change
#include "layered.hpp"                                // layer_h

#ifdef ELECTROSTATICS

namespace Coulomb {

    void pressure_n_coulomb(int &n_coulomb) {
      switch (coulomb.method) {
        case COULOMB_NONE:
          n_coulomb = 0;
          break;
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          n_coulomb = 2;
          break;
        default:
          n_coulomb = 1;
      }
    }

    void pressure_calc_long_range_force(Observable_stat &virials,
                                        Observable_stat &p_tensor) {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
          fprintf(stderr,
                  "WARNING: pressure calculated, but ELC pressure not implemented\n");
          break;
        case COULOMB_P3M_GPU:
          fprintf(
                  stderr,
                  "WARNING: pressure calculated, but GPU P3M pressure not implemented\n");
          break;
        case COULOMB_P3M: {
          p3m_charge_assign();
          virials.coulomb[1] = p3m_calc_kspace_forces(0, 1);
          p3m_charge_assign();
          p3m_calc_kspace_stress(p_tensor.coulomb + 9);
          break;
        }
#endif
        case COULOMB_MMM2D:
          fprintf(
                  stderr,
                  "WARNING: pressure calculated, but MMM2D pressure not implemented\n");
          break;
        case COULOMB_MMM1D:
        case COULOMB_MMM1D_GPU:
          fprintf(
                  stderr,
                  "WARNING: pressure calculated, but MMM1D pressure not implemented\n");
          break;
        default:
          break;
      }
    }

    void nonbonded_coulomb_sanity_checks(int &state) {
      switch (coulomb.method) {
        case COULOMB_MMM1D:
          if (MMM1D_sanity_checks())
            state = 0;
          break;
        case COULOMB_MMM2D:
          if (MMM2D_sanity_checks())
            state = 0;
          break;
#ifdef P3M
        case COULOMB_ELC_P3M:
          if (ELC_sanity_checks())
            state = 0; // fall through
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          if (p3m_sanity_checks())
            state = 0;
          break;
#endif
        default:
          break;
      }
    }

    void nonbonded_electrostatics_cutoff(double &ret) {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
          ret = std::max(elc_params.space_layer, p3m.params.r_cut_iL * box_l[0]);
          break;
        case COULOMB_MMM2D:
          ret = layer_h - skin;
          break;
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          /* do not use precalculated r_cut here, might not be set yet */
          ret = p3m.params.r_cut_iL * box_l[0];
          break;
#endif
        case COULOMB_DH:
          ret = dh_params.r_cut;
          break;
        case COULOMB_RF:
          ret = rf_params.r_cut;
          break;
#ifdef SCAFACOS
        case COULOMB_SCAFACOS:
          ret = Scafacos::get_r_cut();
#endif
        default:
          break;
      }
    }

    void nonbonded_deactivate_coulomb() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          break;
#endif
        case COULOMB_DH:
          dh_params.r_cut = 0.0;
          dh_params.kappa = 0.0;
        case COULOMB_RF:
          rf_params.kappa = 0.0;
          rf_params.epsilon1 = 0.0;
          rf_params.epsilon2 = 0.0;
          rf_params.r_cut = 0.0;
          rf_params.B = 0.0;
        case COULOMB_MMM1D:
          mmm1d_params.maxPWerror = 1e40;
        case COULOMB_MMM2D:
          mmm2d_params.far_cut = 0;
        default:
          break;
      }
    }

    void integrate_coulomb_sanity_check() {
      switch (coulomb.method) {
        case COULOMB_NONE:
          break;
        case COULOMB_DH:
          break;
        case COULOMB_RF:
          break;
#ifdef P3M
        case COULOMB_P3M:
          break;
#endif /*P3M*/
        default: {
          runtimeErrorMsg()
                  << "npt only works with P3M, Debye-Huckel or reaction field";
        }
      }
    }

    void on_observable_calc() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          EVENT_TRACE(
                  fprintf(stderr, "%d: p3m_count_charged_particles\n", this_node));
          p3m_count_charged_particles();
          break;
#endif
        case COULOMB_MAGGS:
          maggs_init();
          break;
        default:
          break;
      }
    }

    void on_coulomb_change() {
      switch (coulomb.method) {
        case COULOMB_DH:
          break;
#ifdef P3M
#ifdef CUDA
        case COULOMB_P3M_GPU:
          p3m_gpu_init(p3m.params.cao, p3m.params.mesh, p3m.params.alpha);
          break;
#endif
        case COULOMB_ELC_P3M:
          ELC_init();
          // fall through
        case COULOMB_P3M:
          p3m_init();
          break;
#endif
        case COULOMB_MMM1D:
          MMM1D_init();
          break;
        case COULOMB_MMM2D:
          MMM2D_init();
          break;
        case COULOMB_MAGGS:
          maggs_init();
          /* Maggs electrostatics needs ghost velocities */
          on_ghost_flags_change();
          break;
        default:
          break;
      }
    }

    void on_resort_particles() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
          ELC_on_resort_particles();
          break;
#endif
        case COULOMB_MMM2D:
          MMM2D_on_resort_particles();
          break;
        default:
          break;
      }
    }

    void on_boxl_change() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
          ELC_init();
          // fall through
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          p3m_scaleby_box_l();
          break;
#endif
        case COULOMB_MMM1D:
          MMM1D_init();
          break;
        case COULOMB_MMM2D:
          MMM2D_init();
          break;
        case COULOMB_MAGGS:
          maggs_init();
          break;
#ifdef SCAFACOS
        case COULOMB_SCAFACOS:
          Scafacos::update_system_params();
          break;
#endif
        default:
          break;
      }
    }

    void init_coulomb() {
      switch (coulomb.method) {
        case COULOMB_DH:
          break;
#ifdef P3M
        case COULOMB_ELC_P3M:
          ELC_init();
          // fall through
        case COULOMB_P3M:
          p3m_init();
          break;
        case COULOMB_P3M_GPU:
          break;
#endif
        case COULOMB_MMM1D:
          MMM1D_init();
          break;
        case COULOMB_MMM2D:
          MMM2D_init();
          break;
        case COULOMB_MAGGS:
          maggs_init();
          /* Maggs electrostatics needs ghost velocities */
          on_ghost_flags_change();
          break;
        default:
          break;
      }
    }

    void calc_long_range_coulomb_force() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
          if (elc_params.dielectric_contrast_on) {
            ELC_P3M_modify_p3m_sums_both();
            ELC_p3m_charge_assign_both();
            ELC_P3M_self_forces();
          } else
            p3m_charge_assign();

          p3m_calc_kspace_forces(1, 0);

          if (elc_params.dielectric_contrast_on)
            ELC_P3M_restore_p3m_sums();

          ELC_add_force();

          break;
#endif
#ifdef CUDA
        case COULOMB_P3M_GPU:
          if (this_node == 0) {
            FORCE_TRACE(printf("Computing GPU P3M forces.\n"));
            p3m_gpu_add_farfield_force();
          }
          /* there is no NPT handling here as long as we cannot compute energies.g
             This is checked in integrator_npt_sanity_checks() when integration
             starts. */
          break;
#endif
#ifdef P3M
        case COULOMB_P3M:
          FORCE_TRACE(printf("%d: Computing P3M forces.\n", this_node));
          p3m_charge_assign();
#ifdef NPT
          if (integ_switch == INTEG_METHOD_NPT_ISO)
            nptiso.p_vir[0] += p3m_calc_kspace_forces(1, 1);
          else
#endif
            p3m_calc_kspace_forces(1, 0);
          break;
#endif
        case COULOMB_MAGGS:
          maggs_calc_forces();
          break;
        case COULOMB_MMM2D:
          MMM2D_add_far_force();
          MMM2D_dielectric_layers_force_contribution();
          break;
#ifdef SCAFACOS
        case COULOMB_SCAFACOS:
          assert(!Scafacos::dipolar());
          Scafacos::add_long_range_force();
          break;
#endif
        default:
          break;
      }
    }

    void energy_calc_long_range_coulomb_energy(Observable_stat &energy) {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_P3M_GPU:
          printf(
                  "long range energy calculation not implemented for GPU P3M\n"); // TODO
          // make
          // right
          break;
        case COULOMB_P3M:
          p3m_charge_assign();
          energy.coulomb[1] = p3m_calc_kspace_forces(0, 1);
          break;
        case COULOMB_ELC_P3M:
          // assign the original charges first
          // they may not have been assigned yet
          p3m_charge_assign();
          if (!elc_params.dielectric_contrast_on)
            energy.coulomb[1] = p3m_calc_kspace_forces(0, 1);
          else {
            energy.coulomb[1] = 0.5 * p3m_calc_kspace_forces(0, 1);
            energy.coulomb[1] += 0.5 * ELC_P3M_dielectric_layers_energy_self();

            //  assign both original and image charges now
            ELC_p3m_charge_assign_both();
            ELC_P3M_modify_p3m_sums_both();

            energy.coulomb[1] += 0.5 * p3m_calc_kspace_forces(0, 1);

            // assign only the image charges now
            ELC_p3m_charge_assign_image();
            ELC_P3M_modify_p3m_sums_image();

            energy.coulomb[1] -= 0.5 * p3m_calc_kspace_forces(0, 1);
          }
          energy.coulomb[2] = ELC_energy();
          break;
#endif
#ifdef SCAFACOS
        case COULOMB_SCAFACOS:
          assert(!Scafacos::dipolar());
          energy.coulomb[1] += Scafacos::long_range_energy();
          break;
#endif
        case COULOMB_MMM2D:
          *energy.coulomb += MMM2D_far_energy();
          *energy.coulomb += MMM2D_dielectric_layers_energy_contribution();
          break;
          /* calculate electric part of energy (only for MAGGS) */
        case COULOMB_MAGGS:
          *energy.coulomb += maggs_electric_energy();
          break;
        default:
          break;
      }
    }

    void energy_n_coulomb(int &n_coulomb) {
      switch (coulomb.method) {
        case COULOMB_NONE:
          n_coulomb = 0;
          break;
        case COULOMB_ELC_P3M:
          n_coulomb = 3;
          break;
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          n_coulomb = 2;
          break;
        case COULOMB_SCAFACOS:
          n_coulomb = 2;
          break;
        default:
          n_coulomb = 1;
      }
    }

    void icc_calc_pair_coulomb_force(Particle *p1, Particle *p2, double *d,
                                     double dist, double dist2, double *force) {
      auto const q1q2 = p1->p.q * p2->p.q;
      if (q1q2 != 0) {
        switch (coulomb.method) {
#ifdef P3M
          case COULOMB_ELC_P3M:
          case COULOMB_P3M_GPU:
          case COULOMB_P3M:
            p3m_add_pair_force(q1q2, d, dist2, dist, force);
            break;
#endif /* P3M */
          case COULOMB_MMM1D:
            add_mmm1d_coulomb_pair_force(q1q2, d, dist2, dist, force);
            break;
          case COULOMB_MMM2D:
            add_mmm2d_coulomb_pair_force(q1q2, d, dist2, dist, force);
            break;
          default:
            break;
        }
      }
    }

    void icc_calc_long_range_force() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M:
          if (elc_params.dielectric_contrast_on) {
            runtimeErrorMsg() << "ICCP3M conflicts with ELC dielectric contrast";
          }
          p3m_charge_assign();
          p3m_calc_kspace_forces(1, 0);
          ELC_add_force();
          break;

#ifdef CUDA
        case COULOMB_P3M_GPU:
          if (this_node == 0) {
            FORCE_TRACE(printf("Computing GPU P3M forces.\n"));
            p3m_gpu_add_farfield_force();
          }
          break;
#endif
        case COULOMB_P3M:
          p3m_charge_assign();
          p3m_calc_kspace_forces(1, 0);
          break;
#endif
        case COULOMB_MMM2D:
          MMM2D_add_far_force();
          MMM2D_dielectric_layers_force_contribution();
          break;
        default:
          break;
      }
    }

    int iccp3m_sanity_check() {
      switch (coulomb.method) {
#ifdef P3M
        case COULOMB_ELC_P3M: {
          if (elc_params.dielectric_contrast_on) {
            runtimeErrorMsg() << "ICCP3M conflicts with ELC dielectric contrast";
            return 1;
          }
          break;
        }
#endif
        case COULOMB_DH: {
          runtimeErrorMsg() << "ICCP3M does not work with Debye-Hueckel.";
          return 1;
        }
        case COULOMB_RF: {
          runtimeErrorMsg() << "ICCP3M does not work with COULOMB_RF.";
          return 1;
        }
        default:
          break;
      }

#ifdef NPT
      if (integ_switch == INTEG_METHOD_NPT_ISO) {
        runtimeErrorMsg() << "ICCP3M does not work in the NPT ensemble";
        return 1;
      }
#endif

      return 0;
    }

    int elc_switch_error() {
#ifdef P3M
      switch (coulomb.method) {
        case COULOMB_P3M_GPU: {
          runtimeErrorMsg()
                  << "ELC tuning failed, ELC is not set up to work with the GPU P3M";
          return ES_ERROR;
        }
        case COULOMB_ELC_P3M:

        case COULOMB_P3M:
          p3m.params.epsilon = P3M_EPSILON_METALLIC;
          coulomb.method = COULOMB_ELC_P3M;
          return ES_OK;
        default:
          break;
      }
#endif
      return ES_ERROR;
    }

    void bcast_coulomb_params() {
      switch (coulomb.method) {
        case COULOMB_NONE:
          // fall through, scafacos has internal parameter propagation
        case COULOMB_SCAFACOS:
          break;
#ifdef P3M
        case COULOMB_ELC_P3M:
          MPI_Bcast(&elc_params, sizeof(ELC_struct), MPI_BYTE, 0, comm_cart);
          // fall through
        case COULOMB_P3M_GPU:
        case COULOMB_P3M:
          MPI_Bcast(&p3m.params, sizeof(p3m_parameter_struct), MPI_BYTE, 0,
                    comm_cart);
          break;
#endif
        case COULOMB_DH:
          MPI_Bcast(&dh_params, sizeof(Debye_hueckel_params), MPI_BYTE, 0, comm_cart);
          break;
        case COULOMB_MMM1D:
        case COULOMB_MMM1D_GPU:
          MPI_Bcast(&mmm1d_params, sizeof(MMM1D_struct), MPI_BYTE, 0, comm_cart);
          break;
        case COULOMB_MMM2D:
          MPI_Bcast(&mmm2d_params, sizeof(MMM2D_struct), MPI_BYTE, 0, comm_cart);
          break;
        case COULOMB_MAGGS:
          MPI_Bcast(&maggs, sizeof(MAGGS_struct), MPI_BYTE, 0, comm_cart);
          break;
        case COULOMB_RF:
        default:
          fprintf(stderr,
                  "%d: INTERNAL ERROR: cannot bcast coulomb params for "
                  "unknown method %d\n",
                  this_node, coulomb.method);
          errexit();
      }
    }

} // namespace Coulomb
#endif // ELECTROSTATICS