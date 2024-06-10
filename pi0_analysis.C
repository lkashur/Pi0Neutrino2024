/**
 * @file pi0_analysis.C
 * @brief ROOT macro to be used with CAFAna executable to run the selection.
 * @author justin.mueller@colostate.edu
*/

#include "pi0_include/analysis.h"
#include "pi0_include/container.h"
#include "pi0_include/csv_maker.h"
#include "sbnana/CAFAna/Core/Binning.h"

using namespace ana;

/**
 * The main function of the selection. Creates a container for the CAFAna
 * Spectrum objects and populates it with a variety of variables that define
 * the selection.
 * @return none.
*/
void pi0_analysis()
{
  ////////////////////////////////////////////
  /// Monte Carlo Analysis
  ////////////////////////////////////////////
  /*
  // Create spectra
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/mc_run2_new_weights/flat/*.flat.root", "spectra_bnb_nu_cosmic_mc_v09_84_00_01.root", -1, 1.9213e19);
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/mc_run2_new_weights/flat/*.flat.root", "spectra_mc_shift_nu2024.root", -1, 1.92082e+19);
  SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/mc_run2_new_weights/flat/*.flat.root", "spectra_mc_nu2024_full.root", -1, 2.68171e+20);
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/bnb_nucosmics_v6.flat.root", "spectra_nucosmics.root", 1.253e19, 2.5e20);
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.flat.root", "spectra_cv.root", -1, 2.5e20);

  // Fill spectra

  // Basic
  spectra.add_spectrum1d("sPi0SubleadingShowerEnergy_All1muNphCut", Binning::Simple(30,0,400), kPi0SubleadingShowerEnergy_All1muNphCut);
  spectra.add_spectrum1d("sPi0LeadingShowerEnergy_All1muNphCut", Binning::Simple(30,0,400), kPi0LeadingShowerEnergy_All1muNphCut);
  spectra.add_spectrum1d("sPi0OpeningAngle_All1muNphCut", Binning::Simple(30, 0, 180), kPi0OpeningAngle_All1muNphCut);
  spectra.add_spectrum1d("sPi0Mass_All1muNphCut", Binning::Simple(50, 0, 500), kPi0Mass_All1muNphCut);
  spectra.add_spectrum1d("sPi0Mass_All1muNphCutCryoE", Binning::Simple(50, 0, 500), kPi0Mass_All1muNphCutCryoE);
  spectra.add_spectrum1d("sPi0Mass_All1muNphCutCryoW", Binning::Simple(50, 0, 500), kPi0Mass_All1muNphCutCryoW);
  spectra.add_spectrum1d("sPi0LeadingShowerConversionDist_All1muNphCut", Binning::Simple(20, 0, 100), kPi0LeadingShowerConversionDist_All1muNphCut);
  spectra.add_spectrum1d("sPi0SubleadingShowerConversionDist_All1muNphCut", Binning::Simple(20, 0, 100), kPi0SubleadingShowerConversionDist_All1muNphCut);
  spectra.add_spectrum1d("sVisibleEnergy_All1muNphCut", Binning::Simple(25, 0, 3000), kVisibleEnergy_All1muNphCut);
  
  // Efficiency
  spectra.add_spectrum2d("sCountTTP_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_NoCut, kCountTTP_NoCut);
  spectra.add_spectrum2d("sCountTTP_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVCut, kCountTTP_FVCut);
  spectra.add_spectrum2d("sCountTTP_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConCut, kCountTTP_FVConCut);
  spectra.add_spectrum2d("sCountTTP_FVConTop1muNphCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_FVConTop1muNphCut, kCountTTP_FVConTop1muNphCut);
  spectra.add_spectrum2d("sCountTTP_All1muNphCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryTTP_All1muNphCut, kCountTTP_All1muNphCut);

  // Purity
  spectra.add_spectrum2d("sCountPTT_NoCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_NoCut, kCountPTT_NoCut);
  spectra.add_spectrum2d("sCountPTT_FVCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVCut, kCountPTT_FVCut);
  spectra.add_spectrum2d("sCountPTT_FVConCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConCut, kCountPTT_FVConCut);
  spectra.add_spectrum2d("sCountPTT_FVConTop1muNphCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_FVConTop1muNphCut, kCountPTT_FVConTop1muNphCut);
  spectra.add_spectrum2d("sCountPTT_All1muNphCut", Binning::Simple(10, 0, 10), Binning::Simple(1, 0, 2), kCategoryPTT_All1muNphCut, kCountPTT_All1muNphCut);

  // Stacked reco quantities
  spectra.add_spectrum2d("sPi0MassPTT_Topology_All1muNphCut", Binning::Simple(10, 0, 10), Binning::Simple(50, 0, 500), kCategoryTopologyPTT_All1muNphCut, kPi0MassPTT_All1muNphCut);
  spectra.add_spectrum2d("sPi0MassPTT_Topology_All1muNphCryoECut", Binning::Simple(10, 0, 10), Binning::Simple(50, 0, 500), kCategoryTopologyPTT_All1muNphCryoECut, kPi0MassPTT_All1muNphCryoECut);
  spectra.add_spectrum2d("sPi0MassPTT_Topology_All1muNphCryoWCut", Binning::Simple(10, 0, 10), Binning::Simple(50, 0, 500), kCategoryTopologyPTT_All1muNphCryoWCut, kPi0MassPTT_All1muNphCryoWCut);
  spectra.add_spectrum2d("sPi0LeadingShowerConversionDistPTT_Topology_All1muNph", Binning::Simple(10, 0, 10), Binning::Simple(20, 0, 100), kCategoryTopologyPTT_All1muNphCut, kPi0LeadingShowerConversionDist_All1muNphCut);
  spectra.add_spectrum2d("sPi0SubleadingShowerConversionDistPTT_Topology_All1muNph", Binning::Simple(10, 0, 10), Binning::Simple(20, 0, 100), kCategoryTopologyPTT_All1muNphCut, kPi0SubleadingShowerConversionDist_All1muNphCut);
  spectra.add_spectrum2d("sVisibleEnergyPTT_Topology_All1muNphCut", Binning::Simple(10, 0, 10), Binning::Simple(25, 0, 3000), kCategoryTopologyPTT_All1muNphCut, kVisibleEnergyPTT_All1muNphCut);

  // Confusion
  spectra.add_spectrum2d("sPrimary_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth, kPrimary);
  spectra.add_spectrum2d("sPID_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth, kPID);
  spectra.add_spectrum2d("sPrimaryPID_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth, kPrimaryPID);
  spectra.add_spectrum2d("sPrimary_Neutrino_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth_Neutrino, kPrimary_Neutrino);
  spectra.add_spectrum2d("sPID_Neutrino_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth_Neutrino, kPID_Neutrino);
  spectra.add_spectrum2d("sPrimaryPID_Neutrino_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth_Neutrino, kPrimaryPID_Neutrino);
  spectra.add_spectrum2d("sPrimary_Cosmic_confusion", Binning::Simple(2,0,2), Binning::Simple(2,0,2), kPrimaryTruth_Cosmic, kPrimary_Cosmic);
  spectra.add_spectrum2d("sPID_Cosmic_confusion", Binning::Simple(5,0,5), Binning::Simple(5,0,5), kPIDTruth_Cosmic, kPID_Cosmic);
  spectra.add_spectrum2d("sPrimaryPID_Cosmic_confusion", Binning::Simple(10,0,10), Binning::Simple(10,0,10), kPrimaryPIDTruth_Cosmic, kPrimaryPID_Cosmic);

  // Output logging
  spectra.add_spectrum1d("sSelected", Binning::Simple(1, 0, 2), kSelected);
  spectra.add_spectrum1d("sSignal", Binning::Simple(1, 0, 2), kSignal);
  //spectra.add_spectrum1d("sPhotonE", Binning::Simple(1, 0, 2), kPhotonE);
  */

  ////////////////////////////////////////////
  /// Data Analysis
  ////////////////////////////////////////////
  
  // Create spectra
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2_new_weights/onbeam/flat/*.flat.root", "spectra_bnb_run2_data_v09_84_00_01.root", -1, 1.9213e+19);
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2_actual_new_weights/onbeam/flat/*.flat.root", "spectra_bnb_run2_data_v09_84_00_01.root", -1, 1.9213e+19);
  //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2_actual_new_weights/onbeam/flat/*.flat.root", "spectra_data_shift_nu2024.root", -1, 1.92082e+19);
  SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/physics_run2_actual_new_weights/offbeam/hdf5/*.flat.root", "spectra_offbeam_shift.root", -1, 266267); // off-beam
  
  // Fill spectra
  
  // Basic
  spectra.add_spectrum1d("sPi0SubleadingShowerEnergy_All1muNphCut", Binning::Simple(30,0,400), kPi0SubleadingShowerEnergy_All1muNphCut);
  spectra.add_spectrum1d("sPi0LeadingShowerEnergy_All1muNphCut", Binning::Simple(30,0,400), kPi0LeadingShowerEnergy_All1muNphCut);               
  spectra.add_spectrum1d("sPi0OpeningAngle_All1muNphCut", Binning::Simple(30, 0, 180), kPi0OpeningAngle_All1muNphCut);                     
  spectra.add_spectrum1d("sPi0Mass_All1muNphCut", Binning::Simple(50, 0, 500), kPi0Mass_All1muNphCut);                                 
  spectra.add_spectrum1d("sPi0Mass_All1muNphCutCryoE", Binning::Simple(50, 0, 500), kPi0Mass_All1muNphCutCryoE);
  spectra.add_spectrum1d("sPi0Mass_All1muNphCutCryoW", Binning::Simple(50, 0, 500), kPi0Mass_All1muNphCutCryoW);
  spectra.add_spectrum1d("sPi0LeadingShowerConversionDist_All1muNphCut", Binning::Simple(20, 0, 100), kPi0LeadingShowerConversionDist_All1muNphCut);
  spectra.add_spectrum1d("sPi0SubleadingShowerConversionDist_All1muNphCut", Binning::Simple(20, 0, 100), kPi0SubleadingShowerConversionDist_All1muNphCut);                          
  spectra.add_spectrum1d("sVisibleEnergy_All1muNphCut", Binning::Simple(25, 0, 3000), kVisibleEnergy_All1muNphCut);
 
  // Output logging
  //spectra.add_spectrum1d("sData", Binning::Simple(1, 0, 2), kData);
  





  // Add systematic info below

    /**
     * 3. BNB neutrino (full flux) + out-of-time cosmics *     Central Value    * (v09_82_02_01).
     * 4. BNB neutrino (full flux) + out-of-time cosmics * Coherent Noise +4.5% * (v09_82_02_01).
     * 5. BNB neutrino (full flux) + out-of-time cosmics *  Elli. Recombination * (v09_82_02_01).
     * 6. BNB neutrino (full flux) + out-of-time cosmics * Untuned Signal Shape * (v09_82_02_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_cv_v2.flat.root", "spectra_cv.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_tpcnoise_coh_p1_v2.flat.root", "spectra_tpcnoise_coh_p1.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_recombination.flat.root", "spectra_recombination.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/ml_hdf5/bnb_nu_sys/systematics_untunedsigshape.flat.root", "spectra_untunedsigshape.root", -1, 2.5e20);
    
    /**
     * 7. BNB neutrino-only (full flux)  *     Central Value    * (v09_82_02_01).
     * 8. BNB neutrino-only (full flux)  * Coherent Noise +4.5% * (v09_82_02_01).
     * 9. BNB neutrino-only (full flux)  *  Elli. Recombination * (v09_82_02_01).
     * 10. BNB neutrino-only (full flux) * Untuned Signal Shape * (v09_82_02_01).
    */
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.flat.root", "spectra_cv.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_cohnoise.flat.root", "spectra_tpcnoise_coh_p1.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_intnoise.flat.root", "spectra_intnoise.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_recombination.flat.root", "spectra_recombination.root", -1, 2.5e20);
    //SpecContainer spectra("/pnfs/icarus/scratch/users/mueller/systematics/sample_sigshape.flat.root", "spectra_untunedsigshape.root", -1, 2.5e20);
    

    spectra.run();
}
