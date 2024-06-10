/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
*/
#ifndef CUTS_H
#define CUTS_H

#include <functional>
#include <vector>
#include <string>
#include <sstream>
#include <numeric>

namespace cuts
{
    /**
     * Apply a cut on whether a match exists.
     * @tparam T the object type (true or reco, interaction or particle).
     * @param obj the object to select on.
     * @return true if the object is matched.
    */
    template<class T>
        bool matched(const T & obj) { return obj.match.size() > 0; }

    /**
     * Apply a cut on the validity of the flash match.
     * @tparam T the type of interaction (true or reco).
     * @param interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
    */
    template<class T>
        bool valid_flashmatch(const T & interaction) { return !std::isnan(interaction.flash_time) && interaction.fmatched == 1; }

    /**
     * Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the count of primaries of each particle type within the interaction.
     */
    template<class T>
      std::vector<uint32_t> count_primaries(const T & interaction)
      {
	std::vector<uint32_t> counts(5, 0);
	for(auto &p : interaction.particles)
          {
	    /*
	    double scalegain;
	    if(p.volume_id == 0 && p.start_point[0] < -210.215)
	      {
		scalegain =  76.02;
	      }
	    else if(p.volume_id == 0 && p.start_point[0] > -210.215)
	      {
		scalegain =  75.29;
	      }
	    else if(p.volume_id == 1 && p.start_point[0] < 210.215)
	      {
		scalegain =  77.36;
	      }
	    else if(p.volume_id == 1 && p.start_point[0] > 210.215)
	      {
		scalegain =  77.10;
	      }
	    */
            if(p.is_primary)
	      {
		double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
		if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			       energy = p.energy_deposit;

		if(p.pid != 0 && p.pid != 1 && p.pid != 2)
		  {
		    counts[p.pid]++;
		  }
		//else if((p.pid == 0 || p.pid == 1) && (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * energy > 40) // MC
		//else if((p.pid == 0 || p.pid == 1) && (1.0238) * (1./77.0777) * (77.0777*(1./0.90)) * (0.87) * (1./0.78) * energy > 40) // shower energy cut (40 MeV) // MC shifted
		else if((p.pid == 0 || p.pid == 1) && (0.9902) * (1./76.44) * (76.44*(1./0.92)) * (0.87) * (1./0.78) * energy > 40) // Data shifted
		  {		    
		    counts[p.pid]++;
		  }
		else if(p.pid == 2 && energy > 143.425) // muon energy cut
		  {
		    counts[p.pid]++;
		  }
     
	      }
          }
	return counts;
      }

    /**
     * Find the topology of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the topology of the interaction as a string (e.g 0ph0e1mu0pi1p).
     */
    template<class T>
        std::string topology(const T & interaction)
        {
            std::vector<uint32_t> counts(count_primaries(interaction));
            std::stringstream ss;
            ss  << counts[0] << "ph"
                << counts[1] << "e"
                << counts[2] << "mu"
                << counts[3] << "pi"
                << counts[4] << "p";
            return ss.str();
        }

    /**
     * Apply no cut (all interactions/particles passed).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to select on.
     * @return true (always).
     */
    template<class T>
        bool no_cut(const T & obj) { return true; }

    /**
     * Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
        bool fiducial_cut(const T & interaction) { return interaction.is_fiducial; }
    
    /**
     * Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
        bool containment_cut(const T & interaction) { return interaction.is_contained; }

    /**
     * Apply an east cryostat cut.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return true if the vertex is in the east cryostat.
     */
    template<class T>
      bool cryostatE_cut(const T & interaction) { return interaction.volume_id == 0; }

    /**
     * Apply a west cryostat cut.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return true if the vertex is in the west cryostat.
     */
    template<class T>
      bool cryostatW_cut(const T & interaction) { return interaction.volume_id == 1; }
    
    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6);
        }

    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB data.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut_data(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= -0.5) && (interaction.flash_time <= 1.4);
        }

    /**
     * Apply a fiducial and containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
     */
    template<class T>
        bool fiducial_containment_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a 1muNph (N>=2) topological cut. The interaction must have a topology
     * matching 1muNph as defined in the conditions of count_primaries().
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction has 1muNph topology.
     */
    template<class T>
      bool topological_1muNph_cut(const T & interaction)
      {
	std::vector<uint32_t> c(count_primaries(interaction));
	return c[0] >= 2 && c[2] == 1;
      }
    
    template<class T>
      bool topological_1mu2ph_cut(const T & interaction)
      {
	std::vector<uint32_t> c(count_primaries(interaction));
	return c[0] == 2 && c[2] == 1;
      }

    /**
     * Apply a 1mu1pi0 topological cut at the truth level.
     * @tparam T the type of interaction (true).
     * @param interaction to select on.
     * @return true if interaction has 1mu1pi0 topology
     */
    template<class T>
      bool topological_1mu1pi0_cut(const T & interaction)
      {
	int primary_muon_count = 0;
	int primary_pi0_count = 0;
	//vector<int> primary_photon_track_ids;
	unordered_map< int, vector<tuple<int, int, float>> > pi0_map;
	for(auto &p : interaction.particles)
	  {
	    // primary muons
	    if(p.pid == 2 && p.is_primary)
	      {
		// Muon energy cut
		if(p.energy_init > 143.425)
		  {
		    primary_muon_count++;
		  }
	      }

	    // primary pi0 photons
	    //if(p.pdg_code == 22 && p.is_primary && p.ancestor_pdg_code == 111 && p.energy_init > 40)
	    //  {
	    //	primary_photon_track_ids.push_back(p.ancestor_track_id);
	    //  }

	    if(p.is_primary && p.ancestor_pdg_code == 111)
	      {
		pi0_map[p.ancestor_track_id].push_back(make_tuple(p.id, p.pdg_code, p.energy_init));
	      }
	      
	  }

	// count duplicate ancestor ids (pairs of photons) to find true pi0s
	//map<int, int> primary_pi0_frequency;
	//for(int i : primary_photon_track_ids)
	//  {
	//    primary_pi0_frequency[i]++;
	//  }
	//primary_pi0_count = count_if(primary_pi0_frequency.begin(), primary_pi0_frequency.end(),
	//[](auto const &d) { return d.second == 2; });  
	
	// Loop over pi0 map
	for (auto const& pi0 : pi0_map)
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
	

	if (primary_muon_count == 1 && primary_pi0_count == 1)
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }
	
      }

    /**
     * Apply a fiducial, containment, and topological (1muNph) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
      bool fiducial_containment_topological_1muNph_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction)\
	&& topological_1muNph_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNph), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction passes the above cuts.
     */
    template<class T>
      bool all_1muNph_cut(const T & interaction) { return topological_1muNph_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_data<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNph), flash time, and east cryostat cut.
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction passes the above cuts.
     */
    template<class T>
      bool all_1muNph_cryoE_cut(const T & interaction) {return all_1muNph_cut(interaction) && cryostatE_cut(interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNph), flash time, and west cryostat cut.
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction passes the above cuts.
     */
    template<class T>
      bool all_1muNph_cryoW_cut(const T & interaction) {return all_1muNph_cut(interaction) && cryostatW_cut(interaction); }

    /**
     * Defined the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
     */
    template<class T>
      bool neutrino(const T & interaction) { return interaction.is_neutrino; }

    /**
     * Define the true cosmic interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a cosmic.
     */
    template<class T>
      bool cosmic(const T & interaction) { return !interaction.is_neutrino; }

    /**
     * Define the true 1mu1pi0 interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muNph neutrino interaction.
     */
    template<class T>
      bool signal_1mu1pi0(const T & interaction) 
      {
	return topological_1mu1pi0_cut<T>(interaction) && neutrino(interaction); 
      }
    
    template<class T>
      bool signal_1mu2ph(const T & interaction)
      {
	return topological_1mu2ph_cut<T>(interaction) && neutrino(interaction);
      }

    /**
     * Define the true "other neutrino" interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
      bool other_nu(const T & interaction) { return !topological_1mu1pi0_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define true photon particle classification.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if particle is a photon.
     */
    template<class T>
      bool photon(const T & particle) { return particle.pid == 0 &&  particle.is_primary && particle.is_contained; }

    /**
     * Define true photon particle classification (matched).
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if particle is a photon (matched).
     */
    template<class T>
      bool matched_photon(const T & particle) { return photon(particle) && matched(particle); }

}
#endif
