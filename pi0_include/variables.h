/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include <algorithm>
#include <TVector3.h>
#include <iostream>

namespace vars
{
  
  /**
   * Variable for image_id (unique identifier for the event).
   * @tparam T the type of object (true or reco, interaction or particle).
   * @param interaction/particle to apply the variable on.
   * @return the image_id of the interaction/particle.
   */
    template<class T>
      double image_id(const T & obj) { return obj.image_id; }

  /**
   * Variable for id (unique identifier for the object).
   * @tparam T the type of object (true or reco, interaction or particle).
   * @param interaction/particle to apply the variable on.
   * @return the id of the interaction/particle.
   */
    template<class T>
      double id(const T & obj) { return obj.id; }

  /**
   * Variable for counting interactions/particles.                                                                                                    
   * @tparam T the type of object (true or reco, interaction or particle).                                                                            
   * @param interaction/particle to apply the variable on.                                                                                            
   * @return 1.0 (always).                                                                                                                            
   */
    template<class T>
      double count(const T & obj) { return 1.0; }

  /**
   * Variable for enumerating interaction categories. This is a basic
   * categorization using only signal, neutrino background, and cosmic
   * background as the three categories.
   * 0: 1mu1pi0 (contained and fiducial)
   * 1: 1mu1pi0 (not contained or fiducial)
   * 2: Other nu interaction
   * 3: cosmic
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the enumerated category of the interaction.
   */
    template<class T>
      double category(const T & interaction)
      {
	double cat(3);
	if(cuts::signal_1mu1pi0(interaction))
	  {
	    if(cuts::fiducial_containment_cut(interaction)) cat = 0;
	    else cat = 1;
	  }
	else if(cuts::other_nu(interaction)) cat = 2;
	return cat;
      }

    /**
     * Variable for particle primary categorizations.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the primary/non-primary designation of the particle.
     */
    template<class T>
      double primary(const T & particle) { return particle.is_primary ? 1 : 0; }

  /**
   * Variable for particle PID.
   * @tparam T the type of particle (true or reco).
   * @param particel to apply the variable on.
   * @return the PID of the particle.
   */
    template<class T>
      double pid(const T & particle) { return particle.pid; }

    /**
     * Variable for particle PID + primary.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID+primary information for the particle.
     */
    template<class T>
      double primary_pid(const T & particle) { return particle.pid + (particle.is_primary ? 5 : 0); }
 
  /**
   * Variable for enumerating interaction categories. This classifies the
   * interactions based on the visible final states.
   * 0: 1mu1pi0, 1: 1muNpi0, 2: 1muCex, 3: NC, 4: Cosmic, 5: Other
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the enumerated category of the interaction.
   */
  template<class T>
    double category_topology(const T & interaction)
    {
      uint16_t cat(4);
      if(interaction.is_neutrino)
	{
	  int primary_muon_count = 0;
	  unordered_map<int, vector<tuple<int, int, float>>> primary_pi0_map;
	  unordered_map<int, vector<tuple<int, int, float>>> cex_pi0_map;
	  int primary_pi0_count = 0;
	  int cex_pi0_count = 0;
	  vector<int> primary_pi0_track_ids;
	  vector<int> cex_pi0_track_ids;
	  for(auto &p : interaction.particles)
	    {
	      // Muon count
	      if(p.pid == 2 && p.is_primary)
		{
		  primary_muon_count++;
		}

	      // Collect primary pi0 photons
	      if(p.is_primary && p.ancestor_pdg_code == 111)
		{
		  primary_pi0_map[p.ancestor_track_id].push_back(make_tuple(p.id, p.pdg_code, p.energy_init));
		}

	      // Collect cex pi0 photons
	      if(!p.is_primary && abs(p.ancestor_pdg_code) == 211)
		{
		  cex_pi0_map[p.ancestor_track_id].push_back(make_tuple(p.id, p.pdg_code, p.energy_init));
		}



	      // Collect primary pi0 photons
	      /*
	      if(p.pdg_code == 22 && p.is_primary && p.ancestor_pdg_code == 111)
		{
		  primary_pi0_track_ids.push_back(p.ancestor_track_id);
		}

	      // Collect cex pi0 photon
	      if(p.pdg_code == 22 && !p.is_primary && abs(p.ancestor_pdg_code) == 211)
		{                                                                                                                                    
		  cex_pi0_track_ids.push_back(p.ancestor_track_id);
		}
	      */

	    }

	  // Number of primary pi0s
	  for (auto const& pi0 : primary_pi0_map)
	    {
	      int num_primary_photon_daughters = 0;
	      for(auto d : pi0.second)
		{
		  if(get<1>(d) == 22 && get<2>(d) > 40) num_primary_photon_daughters++;
		}

	      if(num_primary_photon_daughters == 2)
		{
		  primary_pi0_count++;
		}
	    }
	  
	  // Number of Cex pi0s
	  for (auto const& pi0 : cex_pi0_map)
	    {
	      int num_cex_photon_daughters = 0;
	      for(auto d : pi0.second)
		{
		  if(get<1>(d) == 22) num_cex_photon_daughters++;
		}
	      
	      if(num_cex_photon_daughters == 2)
		{
		  cex_pi0_count++;
		}
	    }


	  // Count duplicate ancestor IDs to find number of primary pi0s
	  /*
	  map<int, int> primary_pi0_frequency;
	  for(int i : primary_pi0_track_ids)
	    {
	      primary_pi0_frequency[i]++;
	    }
	  primary_pi0_count = count_if(primary_pi0_frequency.begin(), primary_pi0_frequency.end(),
				       [](auto const &d) { return d.second == 2; });


	  // Count duplicate ancestor IDs to find number of Cex pi0s
	  map<int, int> cex_pi0_frequency;
	  for(int i : cex_pi0_track_ids)
	    {
	      cex_pi0_frequency[i]++;
	    }
	  cex_pi0_count = count_if(cex_pi0_frequency.begin(), cex_pi0_frequency.end(),
				   [](auto const &d) { return d.second >= 2; });
	  */
	  
	  if(primary_muon_count == 1 && primary_pi0_count == 1 && cex_pi0_count == 0 && interaction.nu_current_type == 0 && interaction.is_contained && interaction.is_fiducial) cat = 0;
	  else if(primary_muon_count == 1 && primary_pi0_count > 1 && cex_pi0_count == 0 && interaction.nu_current_type == 0) cat = 1;
	  else if(primary_muon_count == 1 && cex_pi0_count >= 1 && interaction.nu_current_type == 0) cat = 2;
	  else if(interaction.nu_current_type == 1) cat = 3;
	  else
	    { 
	      cat = 5;
	      //std::cout << "Interaction" << endl;
	      //for(auto &p : interaction.particles)
	      //  {
	      //std::cout << p.is_primary << " " << p.pid << " " <<p.pdg_code << " " << p.ancestor_pdg_code << " " << p.ancestor_track_id << " " << interaction.is_contained << endl;
	      //}
	    }
	}

      return cat;
    }


  /**
   * Variable for total visible energy of interaction.
   * @tparam T the type of interaction (true or reco).
   * @param inteaction to apply the variable on.
   * @return the total visible energy of the interaction.
   */
  template<class T>
    double visible_energy(const T & interaction)
    {
      double energy(0);
      for(const auto & p : interaction.particles)
	{
	  if(p.is_primary)
	    {
	      if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			     {
			       energy += p.energy_deposit;
			     }
	      else
		{
		  if(p.pid < 2) energy += p.calo_ke;
		  else energy += p.csda_ke;
		}
	      if(p.pid == 2) energy += MUON_MASS;
	      else if(p.pid == 3) energy += PION_MASS;
	    }
	}
      return energy;
    }

  /**
   * Variable for energy of the neutrino primary of the interaction.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the neutrino energy.
   */
    template<class T>
      double neutrino_energy(const T & interaction) { return 1000*interaction.nu_energy_init; }

  /**
   * Variable for true particle starting kinetic energy.
   * @tparam T the type of particle (true or reco).
   * @param particle to apply the variable on.
   * @return the starting kinetic kinetic energy of the particle.
   */
    template<class T>
      double ke_init(const T & particle)
      {
	double energy(particle.energy_init);
	switch (particle.pid)
	  {
	  case 1:
	    energy -= ELECTRON_MASS;
	    break;
	  case 2:
	    energy -= MUON_MASS;
	    break;
	  case 3:
	    energy -= PION_MASS;
	    break;
	  case 4:
	    energy -= PROTON_MASS;
	    break;
	  default:
	    break;
	  }
	return energy;
      }

  /**
   * Variable fo particle calo_ke.
   * @tparam T the type of particle (true or reco).
   * @param particle to apply the variable on.
   * @return the calo_ke of the particle.
   */
    template<class T>
      double calo_ke(const T & particle) { return particle.calo_ke; }

  /**
   * Variable for particle calo_ke (photons only).
   * @tparam T the type of particle (true or reco).
   * @param particle to apply the variable on.
   * @return the calo_ke of the particle (if a photon).
   */
    template<class T>
      double calo_ke_photon(const T & particle) { return (cuts::photon(particle)) ? calo_ke(particle) : -1; }

    // temp scaling for incorrect gain in data
    template<class T>
      double tempscale(const T & particle)
      {
	if(particle.volume_id == 0 && particle.start_point[0] < -210.215)
	  {
	    return 76.02;
	  }
	else if(particle.volume_id == 0 && particle.start_point[0] > -210.215)
	  {
	    return 75.29;
	  }
	else if(particle.volume_id == 1 && particle.start_point[0] < 210.215)
	  {
	    return 77.36;
	  }
	else if(particle.volume_id == 1 && particle.start_point[0] > 210.215)
	  {
	    return 77.10;
	  }

      }

  /**
   * Variable for finding true pi0 photons.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return vector true pi0 photon IDs.
   */    
    template<class T>
      vector<int> true_pi0_photon_ids(const T & interaction)
      {
	// Output
	vector<int> pi0_photon_ids;
	unordered_map<int, vector<tuple<int, int, float>>> pi0_map;

	// Signal events
	if(cuts::signal_1mu1pi0(interaction) && cuts::fiducial_containment_cut(interaction))
	  {
	    for(auto & p : interaction.particles)
	      {
		if(p.is_primary && p.ancestor_pdg_code == 111)
		  {
		    pi0_map[p.ancestor_track_id].push_back(make_tuple(p.id, p.pdg_code, p.energy_init));
		  }
	      }

	    // Primary pi0s
	    for (auto const& pi0 : pi0_map)
	      {
		int num_primary_photon_daughters = 0;
		for(auto d : pi0.second)
		  {
		    if(get<1>(d) == 22 && get<2>(d) > 40) num_primary_photon_daughters++;
		  }

		if(num_primary_photon_daughters == 2)
		  {
		    for(auto d : pi0.second)
		      {
			if(get<1>(d) == 22 && get<2>(d) > 40)
			  {
			    pi0_photon_ids.push_back(get<0>(d));
			  }
		      }
		  }
	      }
	  }
	return pi0_photon_ids;
      }
    
    template<class T>
      int true_pi0_leading_photon_id(const T & interaction)
      {

	vector<int> photon_ids = true_pi0_photon_ids(interaction);
	int photon0_id;
	double photon0_energy;
	int photon1_id;
	double photon1_energy;
	// Loop over particles
	for(auto & p : interaction.particles)
	  {
	    if (p.id == photon_ids[0])
	      {
		photon0_id = p.id;
		photon0_energy = p.energy_init;
	      }
	    else if (p.id == photon_ids[1])
	      {
		photon1_id = p.id;
		photon1_energy = p.energy_init;
	      }
	  }

	// Get leading ID
	int leading_id;
	if(photon0_energy > photon1_energy)
	  {
	    leading_id = photon0_id;
	  }
	else
	  {
	    leading_id = photon1_id;
	  }

	return leading_id;
      }

    template<class T>
      int true_pi0_subleading_photon_id(const T & interaction)
      {
	vector<int> photon_ids = true_pi0_photon_ids(interaction);
        int photon0_id;
        double photon0_energy;
        int photon1_id;
        double photon1_energy;

	// Loop over particles
	for(auto & p : interaction.particles)
	  {
	    if (p.id == photon_ids[0])
              {
		photon0_id = p.id;
                photon0_energy = p.energy_init;
              }
            else if (p.id == photon_ids[1])
              {
                photon1_id = p.id;
                photon1_energy = p.energy_init;
              }
	  }
	
	// Get suleading ID
	int subleading_id;
	if(photon0_energy > photon1_energy)
	  {
	    subleading_id = photon1_id;
	  }
	else
	  {
	    subleading_id = photon0_id;
	  }

	return subleading_id;
      }
    
  template<class T>
    double true_pi0_leading_photon_energy(const T & interaction)
    {
      int leading_id = true_pi0_leading_photon_id(interaction);
      double leading_energy;
      for(auto & p : interaction.particles)
        {
          if (p.id == leading_id)
            {
              leading_energy = p.energy_init;
            }
        }
      return leading_energy;
    }

    template<class T>
      double true_pi0_subleading_photon_energy(const T & interaction)
      {
	int subleading_id = true_pi0_subleading_photon_id(interaction);
	double subleading_energy;
	for(auto & p : interaction.particles)
	  {
	    if (p.id == subleading_id)
	      {
		subleading_energy = p.energy_init;
	      }
	  }
	return subleading_energy;
      }

  /**
   * Variable for finding candidate pi0 photons.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return vector of pairs where first element corresponds
   * to shower indices and second element corresponds to
   * angular agreement between showers.
   */
  template<class T>
    vector<pair<pair<int, int>, float>> pi0_photon_candidates(const T & interaction)
    {
      // Output storage
      vector<pair< pair<int, int>, float> > photons;
      for(auto & p : interaction.particles)
	{
	  // primary photons only
	  //if (p.pid != 0 || !p.is_primary || ((1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke) < 40) continue; // MC
	  //if (p.pid != 0 || !p.is_primary || ((1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke) < 40) continue; // MC shifted
	  //double scalegain0 = tempscale(p);
	  if (p.pid != 0 || !p.is_primary || ((0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke) < 40) continue; // Data shifted

	  TVector3 sh0_start(p.start_point[0], p.start_point[1], p.start_point[2]);
	  TVector3 sh0_dir(p.start_dir[0], p.start_dir[1], p.start_dir[2]);

	  // Loop through particles, again
	  for(auto & q : interaction.particles)
	    {
	      // (but don't loop through first particle again)
	      if (p.id == q.id) continue;

	      // Primary photons
	      //if (q.pid != 0 || !q.is_primary || ((1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * q.calo_ke) < 40) continue; // MC
	      //if (q.pid != 0 || !q.is_primary || ((1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * q.calo_ke) < 40) continue; //MC shifted
	      //double scalegain1 = tempscale(q);
	      if (q.pid != 0 || !q.is_primary || ((0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * q.calo_ke) < 40) continue; // Data shifted

	      TVector3 sh1_start(q.start_point[0], q.start_point[1], q.start_point[2]);
	      TVector3 sh1_dir(q.start_dir[0], q.start_dir[1], q.start_dir[2]);

	      // Find "vertex" for each photon pair
	      TVector3 c0;
	      TVector3 c1;
	      TVector3 n = sh0_dir.Cross(sh1_dir);
	      TVector3 n0 = sh0_dir.Cross(n);
	      TVector3 n1 = sh1_dir.Cross(n);
	      float s0 = (sh1_start - sh0_start).Dot(n1) / sh0_dir.Dot(n1);
	      float s1 = (sh0_start - sh1_start).Dot(n0) / sh1_dir.Dot(n0);

	      if (s0 > 0 && s1 > 0)
		{
		  c0 = sh0_start;
		  c1 = sh1_start;
		}
	      else if (s0 > 0 && s1 < 0)
		{
		  c0 = sh0_start;
		  c1 = sh0_start;
		}
	      else if (s0 < 0 && s1 > 0)
		{
		  c0 = sh1_start;
		  c1 = sh1_start;
		}
	      else
		{
		  c0 = sh0_start + s0*sh0_dir;
		  c1 = sh1_start + s1*sh1_dir;
		}

	      float d0 = (sh0_start - c0).Mag();
	      float d1 = (sh1_start - c1).Mag();

	      TVector3 vertex;
	      if (d0 == 0 || d1 == 0)
		{
		  float vertex_x = (c0[0] + c1[0]) / 2;
		  float vertex_y = (c0[1] + c1[1]) / 2;
		  float vertex_z = (c0[2] + c1[2]) / 2;
		  vertex.SetX(vertex_x);
		  vertex.SetY(vertex_y);
		  vertex.SetZ(vertex_z);
		}
	      else
		{
		  float vertex_x = ((c0[0] * d1) + (c1[0] * d0)) / (d1 + d0);
		  float vertex_y = ((c0[1] * d1) + (c1[1] * d0)) / (d1 + d0);
		  float vertex_z = ((c0[2] * d1) + (c1[2] * d0)) / (d1 + d0);
		  vertex.SetX(vertex_x);
		  vertex.SetY(vertex_y);
		  vertex.SetZ(vertex_z);
		}	      
	      
	      // Find mean angular displacement between
	      // <sh_start from vertex> and <sh_dir>
	      float r0 = (sh0_start - vertex).Mag();
	      float r1 = (sh1_start - vertex).Mag();
	      float angle = 0;
	      if (r0 > 0)
		{
		  TVector3 v0 = (sh0_start - vertex).Unit();
		  angle += acos( sh0_dir.Dot(v0) )/2;
		}
	      if (r1 > 0)
		{
		  TVector3 v1 = (sh1_start - vertex).Unit();
		  angle += acos( sh1_dir.Dot(v1) )/2;
		}

	      photons.push_back(make_pair(make_pair(p.id, q.id), angle));
	      sort(photons.begin(), photons.end(), [](const pair<pair<int, int>, float> &a, const pair<pair<int, int>, float> &b)
		   { return a.second < b.second;});	      

	    } // end inner particle loop
	} // end outer particle loop

      return photons;
    }

  /**
   * Variable for finding candidate pi0 photons' angular agreement.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the candidate pi0 photons' angular agreement.
   */
  template<class T>
    double pi0_angular_agreement(const T & interaction)
    {
      // Get photon pairs in interaction
      vector<pair< pair<int, int>, float> > photons = pi0_photon_candidates(interaction);
      double a = photons[0].second;

      return a;
      
    }

  /**
   * Variable for finding candidate muon ID.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the muon ID.
   */
  template<class T>
    int muon_id(const T & interaction)
    {
      
      int muon_id;
      // Loop through particles
      for(auto & p : interaction.particles)
	{
	  if ((p.pid == 2) & (p.is_primary))
	    {
	      muon_id = p.id;
	    }
	}

      return muon_id;
    }



  /**
   * Variable for finding candidate pi0 leading photon ID.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the leading photon ID.
   */
  template<class T>
    int pi0_leading_shower_id(const T & interaction)
    {
      // Get photon pairs in interaction
      vector<pair< pair<int, int>, float> > photons = pi0_photon_candidates(interaction);
      pair<int, int> ph_pair_ids = photons[0].first;
      
      // Loop through particles
      int photon0_id;
      double photon0_energy;
      int photon1_id;
      double photon1_energy;
      for(auto & p : interaction.particles)
	{
	  //double scalegain = tempscale(p);
	  if (p.id == ph_pair_ids.first)
	    {
	      photon0_id = p.id;
	      //photon0_energy = (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC
	      //photon0_energy = (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC shifted
	      photon0_energy = (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke; // Data shifted
	    }
	  else if (p.id == ph_pair_ids.second)
	    {
	      photon1_id = p.id;
	      //photon1_energy = (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC
	      //photon1_energy = (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC shifted
	      photon1_energy = (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke; // Data shifted
	    }
	}
      
      // Get leading ID
      int leading_id;
      if(photon0_energy > photon1_energy)
	{
	  leading_id = photon0_id;
	}
      else
	{
	  leading_id = photon1_id;
	}
      
      return leading_id;
    }
  
  /**
   * Variable for finding candidate pi0 subleading photon ID. 
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the subleading photon ID.
   */
  template<class T>
    int pi0_subleading_shower_id(const T & interaction)
    {
      // Get photon pairs in interaction
      vector<pair< pair<int, int>, float> > photons = pi0_photon_candidates(interaction);
      pair<int, int> ph_pair_ids = photons[0].first;

      // Loop through particles
      int photon0_id;
      double photon0_energy;
      int photon1_id;
      double photon1_energy;
      for(auto & p : interaction.particles)
        {
	  //double scalegain = tempscale(p);
          if (p.id == ph_pair_ids.first)
            {
              photon0_id = p.id;
	      //photon0_energy = (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC
	      //photon0_energy = (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC shifted
	      photon0_energy = (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke; // Data shifted
            }
          else if (p.id == ph_pair_ids.second)
            {
              photon1_id = p.id;
	      //photon1_energy = (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC
              //photon1_energy = (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC shifted
	      photon1_energy = (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke; // Data shifted
            }
        }

      // Get subleading ID
      int subleading_id;
      if(photon0_energy > photon1_energy)
        {
          subleading_id = photon1_id;
        }
      else
        {
          subleading_id = photon0_id;
        }

      return subleading_id;
    }

  /**
   * Variable for finding pi0 leading shower energy.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the leading shower energy in MeV.
   */
  template<class T>
    double pi0_leading_shower_energy(const T & interaction)
    {
      int leading_id = pi0_leading_shower_id(interaction);
      double leading_energy;
      for(auto & p : interaction.particles)
	{
	  //double scalegain = tempscale(p);
	  if (p.id == leading_id)
	    {
	      //leading_energy = (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC
	      //leading_energy = (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC shifted
	      leading_energy = (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke; // Data shifted
	    }
	}
      return leading_energy;
    }
  
  /**
   * Variable for finding pi0 subleading shower energy.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the subleading shower energy in MeV.
   */
  template<class T>
    double pi0_subleading_shower_energy(const T & interaction)
    {
      int subleading_id = pi0_subleading_shower_id(interaction);
      double subleading_energy;
      for(auto & p : interaction.particles)
        {
	  //double scalegain = tempscale(p);
          if (p.id == subleading_id)
            {
	      //subleading_energy = (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC
              //subleading_energy = (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * p.calo_ke; // MC shifted
	      subleading_energy = (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * p.calo_ke; // Data shifted
            }
        }
      return subleading_energy;
    }  

  /**
   * Variable for finding pi0 leading shower start point.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the leading shower start point.
   */
  template<class T> 
    vector<double> pi0_leading_shower_start(const T & interaction)
    {
      int leading_id = pi0_leading_shower_id(interaction);
      vector<double> leading_start(3);
      for(auto & p : interaction.particles)
	{
	  if (p.id == leading_id)
	    {
	      leading_start[0] = p.start_point[0];
	      leading_start[1] = p.start_point[1];
	      leading_start[2] = p.start_point[2];
	    }
	}
      return leading_start;
    }

  /**
   * Variable for finding pi0 subleading shower start point.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the subleading shower start point.
   */
  template<class T>
    vector<double> pi0_subleading_shower_start(const T & interaction)
    {
      int subleading_id = pi0_subleading_shower_id(interaction);
      vector<double> subleading_start(3);
      for(auto & p : interaction.particles)
	{
          if (p.id == subleading_id)
            {
	      subleading_start[0] = p.start_point[0];
	      subleading_start[1] = p.start_point[1];
	      subleading_start[2] = p.start_point[2];
            }
	}
      return subleading_start;
    }

  /**
   * Variable for finding pi0 leading shower start distance to vertex.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the leading shower distance to vertex
   */
  template<class T>
    double pi0_leading_shower_start_to_vertex(const T & interaction)
    {
      // Interaction vertex
      TVector3 int_vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);

      // Showers
      int leading_id = pi0_leading_shower_id(interaction);
      TVector3 leading_dir;
      int subleading_id = pi0_subleading_shower_id(interaction);
      TVector3 subleading_dir;

      TVector3 sh_start;
      for(auto & p : interaction.particles)
	{
	  if(p.id == leading_id)
	    {
	      sh_start.SetX(p.start_point[0]);
	      sh_start.SetY(p.start_point[1]);
	      sh_start.SetZ(p.start_point[2]);
	    }
	}
      
      double d = (sh_start - int_vertex).Mag();
      return d;
    }

  /**
   * Variable for finding pi0 subleading shower start distance to vertex.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the leading shower distance to vertex
   */
    template<class T>
      double pi0_subleading_shower_start_to_vertex(const T & interaction)
      {
	// Interaction vertex
	TVector3 int_vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);

	// Showers                                                                                                                                                                                           
	int leading_id = pi0_leading_shower_id(interaction);
	TVector3 leading_dir;
	int subleading_id = pi0_subleading_shower_id(interaction);
	TVector3 subleading_dir;

	TVector3 sh_start;
	for(auto & p : interaction.particles)
	  {
	    if(p.id == subleading_id)
	      {
		sh_start.SetX(p.start_point[0]);
		sh_start.SetY(p.start_point[1]);
		sh_start.SetZ(p.start_point[2]);
	      }
	  }

	double d = (sh_start - int_vertex).Mag();
	return d;
      }

  /**
   * Variable for finding pi0 opening angle.
   * @tparam T the type of interaction (true or reco).
   * @param interaction to apply the variable on.
   * @return the pi0 opening angle.
   */
  template<class T>
    double pi0_opening_angle(const T & interaction)
    {
      // Interaction vertex
      TVector3 int_vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);
      
      int leading_id = pi0_leading_shower_id(interaction);
      TVector3 leading_dir;
      int subleading_id = pi0_subleading_shower_id(interaction);
      TVector3 subleading_dir;
      for(auto & p : interaction.particles)
	{
	  if (p.id == leading_id)
	    {
              leading_dir.SetX(p.start_point[0] - int_vertex[0]);
              leading_dir.SetY(p.start_point[1] - int_vertex[1]);
              leading_dir.SetZ(p.start_point[2] - int_vertex[2]);
              leading_dir = leading_dir.Unit();	      
	    }
	  else if (p.id == subleading_id)
	    {
	      subleading_dir.SetX(p.start_point[0] - int_vertex[0]);
              subleading_dir.SetY(p.start_point[1] - int_vertex[1]);
              subleading_dir.SetZ(p.start_point[2] - int_vertex[2]);
              subleading_dir = subleading_dir.Unit();
	    }
	  else
	    {
	      continue;
	    }
	}

      double cos_opening_angle = leading_dir.Dot(subleading_dir);
      double opening_angle = acos(cos_opening_angle) * (180 / 3.14159);
      return opening_angle;
    }

  template<class T>
    double pi0_mass(const T & interaction)
    {
      double ph0_energy = pi0_subleading_shower_energy(interaction);
      double ph1_energy = pi0_leading_shower_energy(interaction);
      double theta_d = pi0_opening_angle(interaction);
      double theta_r = theta_d * (3.14159 / 180);
      double mass = sqrt(2*ph0_energy*ph1_energy*(1-cos(theta_r)));
      return mass;
    }    

}

#endif
