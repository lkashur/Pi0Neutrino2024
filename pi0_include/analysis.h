/**
 * @file analysis.h
 * @brief Header file for definitions of complete variables with relevant cuts.
 * @author justin.mueller@colostate.edu
*/

#include "definitions.h"
#include "cuts.h"
#include "variables.h"

using namespace ana;

namespace ana
{
    // "Simple" variables that attach to a single interaction attribute.
    VARDLP_RECO(kPi0Mass_All1muNphCut,vars::pi0_mass,cuts::all_1muNph_cut);
    VARDLP_RECO(kPi0Mass_All1muNphCutCryoE,vars::pi0_mass,cuts::all_1muNph_cryoE_cut);
    VARDLP_RECO(kPi0Mass_All1muNphCutCryoW,vars::pi0_mass,cuts::all_1muNph_cryoW_cut);
    VARDLP_RECO(kPi0SubleadingShowerEnergy_All1muNphCut,vars::pi0_subleading_shower_energy,cuts::all_1muNph_cut);
    VARDLP_RECO(kPi0LeadingShowerEnergy_All1muNphCut,vars::pi0_leading_shower_energy,cuts::all_1muNph_cut);
    VARDLP_RECO(kPi0LeadingShowerConversionDist_All1muNphCut,vars::pi0_leading_shower_start_to_vertex,cuts::all_1muNph_cut);
    VARDLP_RECO(kPi0SubleadingShowerConversionDist_All1muNphCut,vars::pi0_subleading_shower_start_to_vertex,cuts::all_1muNph_cut);
    VARDLP_RECO(kPi0OpeningAngle_All1muNphCut,vars::pi0_opening_angle,cuts::all_1muNph_cut);
    VARDLP_RECO(kPi0AngularAgreement_All1muNphCut,vars::pi0_opening_angle,cuts::all_1muNph_cut);
    VARDLP_RECO(kVisibleEnergy_All1muNphCut,vars::visible_energy,cuts::all_1muNph_cut);

    // Define "variables" for binning interactions by some categorical class.
    DEFINECAT();

    // Define variables that are broadcasted across each selection cut.
    TCATVAR(kCountTTP,count);
    TCATVAR(kVisibleEnergyTTP,visible_energy);
    RCATVAR(kVisibleEnergyPTT,visible_energy);
    RCATVAR(kCountPTT,count);
    RCATVAR(kPi0LeadingShowerConversionDistPTT, pi0_leading_shower_start_to_vertex);
    RCATVAR(kPi0SubleadingShowerConversionDistPTT, pi0_subleading_shower_start_to_vertex);
    RCATVAR(kPi0MassPTT, pi0_mass);
    
    // Variables for confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth,vars::primary,cuts::no_cut,cuts::matched);
    PVARDLP_TRUE(kPIDTruth,vars::pid,cuts::no_cut,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth,vars::primary_pid,cuts::no_cut,cuts::matched);
    PVAR_TTP(kPrimary,vars::primary,cuts::no_cut,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID,vars::pid,cuts::no_cut,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID,vars::primary_pid,cuts::no_cut,cuts::no_cut,cuts::no_cut);

    // Variables for neutrino-only confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth_Neutrino,vars::primary,cuts::neutrino,cuts::matched);
    PVARDLP_TRUE(kPIDTruth_Neutrino,vars::pid,cuts::neutrino,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth_Neutrino,vars::primary_pid,cuts::neutrino,cuts::matched);
    PVAR_TTP(kPrimary_Neutrino,vars::primary,cuts::neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID_Neutrino,vars::pid,cuts::neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Neutrino,vars::primary_pid,cuts::neutrino,cuts::no_cut,cuts::no_cut);

    // Variables for cosmic-only confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth_Cosmic,vars::primary,cuts::cosmic,cuts::matched);
    PVARDLP_TRUE(kPIDTruth_Cosmic,vars::pid,cuts::cosmic,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth_Cosmic,vars::primary_pid,cuts::cosmic,cuts::matched);
    PVAR_TTP(kPrimary_Cosmic,vars::primary,cuts::cosmic,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID_Cosmic,vars::pid,cuts::cosmic,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Cosmic,vars::primary_pid,cuts::cosmic,cuts::no_cut,cuts::no_cut);    

    // Variables for 2D "true vs. reco" style plots.
    PVARDLP_TRUE(kCaloKETruth_photon,vars::ke_init,cuts::neutrino,cuts::matched_photon);
    PVAR_TTP(kCaloKE_photon,vars::calo_ke,cuts::neutrino,cuts::photon,cuts::no_cut);

}
