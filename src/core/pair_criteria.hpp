#ifndef PAIR_CRITERIA_HPP
#define PAIR_CRITERIA_HPP

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "energy_inline.hpp"




/** @brief Criterion which provides a true/false for a pair of particles */
class PairCriterion {
  public: 
    virtual bool decide(const Particle& p1, const Particle& p2) =0;
};

/** @brief True if two particles are closer than a cut off distance, respecting minimum image convention */
class DistanceCriterion : public PairCriterion {
  public: 
    DistanceCriterion(double _cut_off) {
      m_cut_off=_cut_off;
    }
    virtual bool decide(const Particle& p1, const Particle& p2) {
      double vec21[3];
      get_mi_vector(vec21,p1.r.p, p2.r.p); 
      return (sqrt(sqrlen(vec21)<= m_cut_off));
    };
    double get_cut_off() {
      return m_cut_off;
    }
    void set_cut_off(double c){
      m_cut_off =c;
    }
    private:
      double m_cut_off;
};


/** True if the short range energy is largern than a cut_off */
class EnergyCriterion : public PairCriterion {
  public: 
    EnergyCriterion(double _cut_off) {
      m_cut_off=_cut_off;
    };
    virtual bool decide(const Particle& p1, const Particle& p2)  {
      double vec21[3];
      const double dist_betw_part =sqrt(distance2vec(p1.r.p, p2.r.p, vec21));
      IA_parameters *ia_params = get_ia_param(p1.p.type, p2.p.type);
      return (calc_non_bonded_pair_energy(const_cast<Particle*>(&p1), const_cast<Particle*>(&p2), ia_params,
                       vec21, dist_betw_part,dist_betw_part*dist_betw_part)) >= m_cut_off;
    };
    double get_cut_off() {
      return m_cut_off;
    }
    void set_cut_off(double c){
      m_cut_off =c;
    }
    private:
      double m_cut_off;
};

/** True if a bond of given type exists between the two particles */
class BondCriterion : public PairCriterion {
  public: 
    BondCriterion(int _bond_type) {
       m_bond_type=_bond_type;
    };
    virtual bool decide(const Particle& p1, const Particle& p2) {
      return bond_exists(&p1,&p2,m_bond_type) || bond_exists(&p2,&p1,m_bond_type);
    };
    int get_bond_type() {
      return m_bond_type;
    };
    void set_bond_type(int t){
      m_bond_type =t;
    }
    private:
      int m_bond_type;
};


#endif

