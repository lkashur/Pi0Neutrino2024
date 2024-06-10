/**
 * @file csv_maker.h
 * @brief Header file defining a dummy SpillMultiVar for dumping particle info.
 * @author justin.mueller@colostate.edu
*/
#ifndef CSV_MAKER_H
#define CSV_MAKER_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "cuts.h"
#include "variables.h"

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

std::ofstream pi0_selection_output("pi0_selection_output.log");
std::ofstream pi0_signal_output("pi0_signal_output.log");
std::ofstream gamma_output("gamma_output.log");

#define GUARD(VAL) std::isinf(VAL) ? -9999 : VAL
#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

/**
 * Writes information about selected interactions to output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the current spill.
 * @return a vector with a single dummy entry.
 */
// PURITY
const SpillMultiVar kSelected([](const caf::SRSpillProxy* sr)
{
  // Loop over reco interactions
  for(auto const & i : sr->dlp)
    {
      // Selected pi0 interactions
      if(cuts::all_1muNph_cut(i))
	{
	  if(cuts::matched(i))
	    {
	      // true matched interaction
	      const auto & j = sr->dlp_true[i.match[0]];
	   
	      pi0_selection_output  //<< CSV((string) sr->hdr.sourceName)
				    << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
				    << CSV(sr->hdr.evt) << CSV(j.nu_id)
				    << CSV(vars::image_id(j)) << CSV(vars::id(j))
				    << CSV(vars::category(j))
				    << CSV(vars::category_topology(j))
				    << CSV(j.vertex[0])
		//<< CSV(cuts::topology(j))
				    << CSV(vars::pi0_subleading_shower_energy(i))
				    << CSV(vars::pi0_leading_shower_energy(i))
				    << CSV(vars::pi0_mass(i))
				    << CSV(vars::visible_energy(i))
				    << CSV(i.vertex[0])
				    << std::endl;
	    } // end matched pi0 loop
	  else
	    {
	      pi0_selection_output  //<< CSV((string) sr->hdr.sourceName)
                                    << CSV(sr->hdr.run) << CSV(sr->hdr.subrun)
                                    << CSV(sr->hdr.evt) << CSV("None")
                                    << CSV("None") << CSV("None")
                                    << CSV("None")
                                    << CSV("None")
				    << CSV("None")
		//<< CSV("None")
				    << CSV(vars::pi0_subleading_shower_energy(i))
                                    << CSV(vars::pi0_leading_shower_energy(i))
                                    << CSV(vars::pi0_mass(i))
				    << CSV(vars::visible_energy(i))
				    << CSV(i.vertex[0])
                                    << std::endl;
	    } // end unmatched pi0 loop
	} //end pi0 interaction loop
    } // end interaction loop

  return std::vector<double>{1};

});

/**
 * Writes information about signal events to output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the current spill.
 * @return a vector with a single dummy entry. 
 */
const SpillMultiVar kSignal([](const caf::SRSpillProxy* sr)
{
  // Loop over true interactions
  for(auto const & ti : sr->dlp_true)
    {
      // 1mu1pi0
      if(cuts::signal_1mu1pi0(ti) && cuts::fiducial_containment_cut(ti))
	{
	  // matched to selected reco interaction
	  if(cuts::matched(ti) && cuts::all_1muNph_cut(sr->dlp[ti.match[0]]))
	    {
	      pi0_signal_output << CSV(vars::image_id(ti))
				<< CSV(vars::id(ti))
				<< CSV(vars::neutrino_energy(ti))
				<< CSV(vars::true_pi0_leading_photon_energy(ti))
				<< CSV(vars::true_pi0_subleading_photon_energy(ti))
				<< CSV("MATCH_SELECTED")
				<< endl;
	    }
	  // matched to other reco interaction
	  else if(cuts::matched(ti) && !cuts::all_1muNph_cut(sr->dlp[ti.match[0]]))
	    {
	      pi0_signal_output << CSV(vars::image_id(ti))
				<< CSV(vars::id(ti))
				<< CSV(vars::neutrino_energy(ti))
				<< CSV(vars::true_pi0_leading_photon_energy(ti))
				<< CSV(vars::true_pi0_subleading_photon_energy(ti))
				<< CSV("MATCH_NOT_SELECTED")
				<< endl;
	    }
	  // not matched to reco interaction
	  else if(!cuts::matched(ti))
	    {
	      pi0_signal_output << CSV(vars::image_id(ti))
				<< CSV(vars::id(ti))
				<< CSV(vars::neutrino_energy(ti))
				<< CSV(vars::true_pi0_leading_photon_energy(ti))
				<< CSV(vars::true_pi0_subleading_photon_energy(ti))
                                << CSV("NO_MATCH")
                                << endl;
	    }
	}
      
    } // end true interaction loop

  return std::vector<double>{1};

}
);


/**
 * Writes information about true, matched photon energy to output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the current spill.
 * @return a vector with a single dummy entry.
 */
const SpillMultiVar kPhotonE([](const caf::SRSpillProxy* sr)
{
 
  // Loop over true interactions                                                                                                             
  for(auto const & ti : sr->dlp_true)
    {

      // Matched reco interaction                                                                                                         
      if(!cuts::matched(ti)) {continue;}
      const auto & ri = sr->dlp[ti.match[0]];

      // Loop over true particles                                                                                                         
      for(auto & tp : ti.particles)
	{
	  // contained photons
	  if(cuts::matched_photon(tp))
	    {
	      // Matched                                                                                                 
	      if(!tp.matched) {continue;}
	      float true_energy_init = tp.energy_init;
	      float true_energy_dep = tp.energy_deposit;
	      float true_calo_ke = 0.87 * tp.calo_ke;
	      float reco_calo_ke;

	      // Loop over reco particles to find match
	      for(auto & rp : ri.particles)
		{
		  if(cuts::matched_photon(rp))
		    {
		      if(rp.id == tp.match[0])
			{
			  reco_calo_ke = (1./77.0777) * (77.0777*(1./0.89)) * 0.87 * rp.calo_ke;
			  gamma_output << CSV(true_energy_init) << CSV(true_energy_dep) << CSV(true_calo_ke) << CSV(reco_calo_ke) << std::endl;
			}
		    }
		}

	    } // matched particle loop
	} // true particle loop
    } // true interaction loop

  return std::vector<double>{1};
});

#endif
