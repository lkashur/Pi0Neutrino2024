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

std::ofstream output("output.log");

#define GUARD(VAL) std::isinf(VAL) ? -9999 : VAL
#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

/**
 * Writes information about true, matched photons to output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the current spill.
 * @return a vector with a single dummy entry.
 */
void kGamma([](const caf::SRSpillProxy* sr)
{
  // Loop over true interactions
  for(auto const & ti : sr->dlp_true)
    {

      // Matched reco interaction
      const auto & ri = sr->dlp[ti.match[0]];

      // Loop over true particles
      for(auto const & tp : ti.particles)
	{
	  // contained photons
	  if(cuts::matched_photon(tp))
	    {
	      // Matched
	      if(tp.matched == False) {continue;}
	      float true_energy_init = tp.energy_init;
	      float true_energy_dep = tp.energy_deposit;
	      float true_calo_ke = 0.87 * tp.calo_ke;
	      float reco_calo_ke;

	      // Loop over reco particles to find match
	      for(auto const & rp : ri.particles)
		{
		  if(cuts::matched_photon(rp))
		    {
		      if(rp.id == tp.match[0])
			{
			  reco_calo_ke = 0.87 * rp.calo_ke;
			  output << CSV(true_energy_init) << CSV(true_energy_dep) << CSV(true_calo_ke) << CSV(reco_calo_ke);
			}
		    }
		}
	      
	    } // matched particle loop
	  
	} // true particle loop
      
    } // true interaction loop

  return std::vector<double>{1};
 });


#endif
