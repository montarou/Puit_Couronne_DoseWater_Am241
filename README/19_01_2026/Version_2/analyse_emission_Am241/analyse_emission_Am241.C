///////////////////////////////////////////////////////////////////////////////
/// @file analyse_emission_Am241.C
/// @brief Analyse statistique de l'émission de la source Am-241 (Geant4)
///
/// Ce script se concentre sur:
/// 1. Distribution du nombre de gammas par événement
/// 2. Distribution du nombre de particules par raie
/// 3. Nombre de gammas primaires atteignant l'eau par événement
/// 4. Histogramme 2D: raie vs multiplicité par événement
/// 5. Gammas ayant interagi dans l'air avant d'atteindre l'eau
///
/// Utilisation:
///   root -l 'analyse_emission_Am241.C("output_100000.root")'
///
/// Auteur: Script pour analyse de simulation Geant4 Am-241
/// Date: Janvier 2026
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TPie.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

// ═══════════════════════════════════════════════════════════════════════════
// CONSTANTES DE LA SIMULATION Am-241
// ═══════════════════════════════════════════════════════════════════════════

const Int_t kNbGammaLines = 12;

const Double_t kGammaLineEnergies[kNbGammaLines] = {
    11.89, 13.9, 17.0, 20.8, 26.3446, 33.1963,
    43.420, 55.56, 59.5409, 98.97, 102.98, 125.30
};

const char* kGammaLineNames[kNbGammaLines] = {
    "X_Ll (11.9 keV)", "X_La (13.9 keV)", "X_Lb (17.0 keV)", "X_Lg (20.8 keV)",
    "g_26keV", "g_33keV", "g_43keV", "g_56keV",
    "g_59.5keV", "g_99keV", "g_103keV", "g_125keV"
};

const char* kGammaLineNamesShort[kNbGammaLines] = {
    "X_Ll", "X_La", "X_Lb", "X_Lg", "g26", "g33",
    "g43", "g56", "g59.5", "g99", "g103", "g125"
};

// Intensités théoriques (%) - Source LNHB
const Double_t kGammaLineIntensities[kNbGammaLines] = {
    1.0, 13.0, 18.5, 5.16, 2.31, 0.1215, 0.0669, 0.0181, 35.92, 0.0203, 0.0195, 0.0041
};

// Moyenne théorique de gammas par désintégration
const Double_t kMeanGammaTheory = 0.7631;


// ═══════════════════════════════════════════════════════════════════════════
// FONCTION PRINCIPALE
// ═══════════════════════════════════════════════════════════════════════════

void analyse_emission_Am241(const char* filename = "output_100000.root") {
    
    // Configuration du style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPalette(kViridis);
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gStyle->SetHistLineWidth(2);
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  ANALYSE STATISTIQUE DE L'ÉMISSION - SOURCE Am-241                   ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Fichier: " << std::setw(58) << std::left << filename << "║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════╝\n\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // OUVERTURE DU FICHIER ROOT
    // ───────────────────────────────────────────────────────────────────────
    
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERREUR: Impossible d'ouvrir le fichier " << filename << std::endl;
        return;
    }
    
    // Charger les arbres
    TTree* eventTree = (TTree*)file->Get("EventData");
    TTree* gammaTree = (TTree*)file->Get("GammaData");
    
    if (!eventTree || !gammaTree) {
        std::cerr << "ERREUR: Arbres EventData ou GammaData non trouvés!" << std::endl;
        return;
    }
    
    Long64_t nEvents = eventTree->GetEntries();
    Long64_t nGammas = gammaTree->GetEntries();
    
    std::cout << ">>> Données chargées: " << nEvents << " événements, " 
              << nGammas << " gammas primaires\n\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // VARIABLES POUR LA LECTURE DES ARBRES
    // ═══════════════════════════════════════════════════════════════════════
    
    // EventData
    Int_t evt_EventID, evt_NGammaEmitted, evt_NGammaWater;
    eventTree->SetBranchAddress("EventID", &evt_EventID);
    eventTree->SetBranchAddress("NGammaEmitted", &evt_NGammaEmitted);
    eventTree->SetBranchAddress("NGammaWater", &evt_NGammaWater);
    
    // GammaData
    Int_t gam_EventID, gam_LineID, gam_ReachedWater, gam_Absorbed;
    Double_t gam_Energy;
    gammaTree->SetBranchAddress("EventID", &gam_EventID);
    gammaTree->SetBranchAddress("Energy", &gam_Energy);
    gammaTree->SetBranchAddress("LineID", &gam_LineID);
    gammaTree->SetBranchAddress("ReachedWater", &gam_ReachedWater);
    gammaTree->SetBranchAddress("Absorbed", &gam_Absorbed);
    
    // ═══════════════════════════════════════════════════════════════════════
    // ANALYSE 1: Distribution du nombre de gammas par événement
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  1. DISTRIBUTION DU NOMBRE DE GAMMAS PAR ÉVÉNEMENT\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    // Histogramme de multiplicité
    TH1D* hNGammaEmitted = new TH1D("hNGammaEmitted", 
        "Nombre de gammas emis par evenement;N_{#gamma} emis;Evenements", 
        15, -0.5, 14.5);
    hNGammaEmitted->SetLineColor(kBlue+1);
    hNGammaEmitted->SetFillColor(kCyan);
    hNGammaEmitted->SetFillStyle(1001);
    
    TH1D* hNGammaWater = new TH1D("hNGammaWater",
        "Nombre de gammas atteignant l'eau;N_{#gamma} #rightarrow eau;Evenements",
        15, -0.5, 14.5);
    hNGammaWater->SetLineColor(kGreen+2);
    hNGammaWater->SetFillColor(kGreen);
    hNGammaWater->SetFillStyle(1001);
    
    // Compteurs par multiplicité
    std::map<Int_t, Long64_t> multCounts;
    std::map<Int_t, Long64_t> multWaterCounts;
    std::map<Int_t, Double_t> sumWaterGivenEmitted; // Somme des gammas→eau pour chaque multiplicité émise
    
    Double_t sumGamma = 0., sumGamma2 = 0.;
    Double_t sumWater = 0., sumWater2 = 0.;
    Int_t maxGamma = 0;
    
    for (Long64_t i = 0; i < nEvents; ++i) {
        eventTree->GetEntry(i);
        
        hNGammaEmitted->Fill(evt_NGammaEmitted);
        hNGammaWater->Fill(evt_NGammaWater);
        
        multCounts[evt_NGammaEmitted]++;
        multWaterCounts[evt_NGammaWater]++;
        sumWaterGivenEmitted[evt_NGammaEmitted] += evt_NGammaWater;
        
        sumGamma += evt_NGammaEmitted;
        sumGamma2 += evt_NGammaEmitted * evt_NGammaEmitted;
        sumWater += evt_NGammaWater;
        sumWater2 += evt_NGammaWater * evt_NGammaWater;
        
        if (evt_NGammaEmitted > maxGamma) maxGamma = evt_NGammaEmitted;
    }
    
    Double_t meanGamma = sumGamma / nEvents;
    Double_t stdGamma = TMath::Sqrt(sumGamma2 / nEvents - meanGamma * meanGamma);
    Double_t meanWater = sumWater / nEvents;
    Double_t stdWater = TMath::Sqrt(sumWater2 / nEvents - meanWater * meanWater);
    
    // Affichage tableau
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ N_gamma │  Nombre d'événements  │  Fraction (%)  │  Cumulé (%)     │\n";
    std::cout << "├─────────┼───────────────────────┼────────────────┼─────────────────┤\n";
    
    Double_t cumul = 0.;
    for (auto& p : multCounts) {
        Double_t frac = 100. * p.second / nEvents;
        cumul += frac;
        std::cout << "│    " << std::setw(2) << p.first << "   │       " 
                  << std::setw(8) << p.second << "        │     " 
                  << std::setw(6) << std::fixed << std::setprecision(2) << frac << "     │     "
                  << std::setw(6) << cumul << "      │\n";
    }
    std::cout << "└─────────┴───────────────────────┴────────────────┴─────────────────┘\n\n";
    
    std::cout << "  Moyenne simulée:   " << std::setprecision(4) << meanGamma << " gamma/événement\n";
    std::cout << "  Moyenne théorique: " << kMeanGammaTheory << " gamma/événement\n";
    std::cout << "  Écart relatif:     " << std::setprecision(2) 
              << 100.*(meanGamma - kMeanGammaTheory)/kMeanGammaTheory << "%\n";
    std::cout << "  Écart-type:        " << std::setprecision(4) << stdGamma << "\n";
    std::cout << "  Maximum observé:   " << maxGamma << " gammas\n\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // ANALYSE 2: Distribution par raie
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  2. DISTRIBUTION DU NOMBRE DE GAMMAS PAR RAIE\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    // Compteurs par raie
    Long64_t lineCounts[kNbGammaLines] = {0};
    Long64_t lineReachedWater[kNbGammaLines] = {0};
    Long64_t lineNotReachedWater[kNbGammaLines] = {0};
    
    // Matrice raie × multiplicité (pour analyse 4)
    // On stocke dans un map: (eventID) -> (NGammaEmitted, liste des LineID)
    std::map<Int_t, Int_t> eventMultiplicity;
    std::map<Int_t, std::vector<Int_t>> eventLineIDs;
    
    // Première passe sur EventData pour récupérer les multiplicités
    for (Long64_t i = 0; i < nEvents; ++i) {
        eventTree->GetEntry(i);
        eventMultiplicity[evt_EventID] = evt_NGammaEmitted;
    }
    
    // Passe sur GammaData
    for (Long64_t i = 0; i < nGammas; ++i) {
        gammaTree->GetEntry(i);
        
        if (gam_LineID >= 0 && gam_LineID < kNbGammaLines) {
            lineCounts[gam_LineID]++;
            
            if (gam_ReachedWater == 1) {
                lineReachedWater[gam_LineID]++;
            } else {
                lineNotReachedWater[gam_LineID]++;
            }
            
            // Stocker pour la matrice 2D
            eventLineIDs[gam_EventID].push_back(gam_LineID);
        }
    }
    
    Long64_t totalGammas = 0;
    for (Int_t l = 0; l < kNbGammaLines; ++l) totalGammas += lineCounts[l];
    
    Double_t sumTheoIntensity = 0.;
    for (Int_t l = 0; l < kNbGammaLines; ++l) sumTheoIntensity += kGammaLineIntensities[l];
    
    std::cout << "┌──────┬─────────────────────────┬────────────┬─────────────┬────────────┬───────────┐\n";
    std::cout << "│  ID  │ Raie                    │ E (keV)    │ N_émis      │ Simu (%)   │ Théo (%)  │\n";
    std::cout << "├──────┼─────────────────────────┼────────────┼─────────────┼────────────┼───────────┤\n";
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        Double_t simuPct = totalGammas > 0 ? 100. * lineCounts[l] / totalGammas : 0.;
        Double_t theoPct = 100. * kGammaLineIntensities[l] / sumTheoIntensity;
        
        std::cout << "│  " << std::setw(2) << l << "  │ " << std::setw(23) << std::left << kGammaLineNamesShort[l] << std::right
                  << " │ " << std::setw(8) << std::setprecision(2) << std::fixed << kGammaLineEnergies[l] << "   │ "
                  << std::setw(9) << lineCounts[l] << "   │ "
                  << std::setw(8) << std::setprecision(3) << simuPct << "   │ "
                  << std::setw(7) << theoPct << "   │\n";
    }
    
    std::cout << "├──────┼─────────────────────────┼────────────┼─────────────┼────────────┼───────────┤\n";
    std::cout << "│ TOT  │                         │            │ " << std::setw(9) << totalGammas 
              << "   │  100.000   │ 100.000   │\n";
    std::cout << "└──────┴─────────────────────────┴────────────┴─────────────┴────────────┴───────────┘\n\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // ANALYSE 3: Gammas atteignant l'eau
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  3. GAMMAS PRIMAIRES ATTEIGNANT L'EAU PAR ÉVÉNEMENT\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ N_γ→eau │  Nombre d'événements  │  Fraction (%)  │  Cumulé (%)     │\n";
    std::cout << "├─────────┼───────────────────────┼────────────────┼─────────────────┤\n";
    
    cumul = 0.;
    for (auto& p : multWaterCounts) {
        Double_t frac = 100. * p.second / nEvents;
        cumul += frac;
        std::cout << "│    " << std::setw(2) << p.first << "   │       " 
                  << std::setw(8) << p.second << "        │     " 
                  << std::setw(6) << std::fixed << std::setprecision(2) << frac << "     │     "
                  << std::setw(6) << cumul << "      │\n";
    }
    std::cout << "└─────────┴───────────────────────┴────────────────┴─────────────────┘\n\n";
    
    std::cout << "  Moyenne:           " << std::setprecision(4) << meanWater << " gamma/événement atteignant l'eau\n";
    std::cout << "  Écart-type:        " << stdWater << "\n";
    std::cout << "  Taux de passage:   " << std::setprecision(2) << 100.*meanWater/meanGamma << "% des gammas émis atteignent l'eau\n\n";
    
    // Tableau comparatif par multiplicité
    std::cout << "  Comparaison par multiplicité:\n";
    std::cout << "  ┌─────────┬────────────┬────────────┬────────────┐\n";
    std::cout << "  │ N_émis  │ N_evt      │ Moy→eau    │ Taux (%)   │\n";
    std::cout << "  ├─────────┼────────────┼────────────┼────────────┤\n";
    
    for (auto& p : multCounts) {
        Int_t nEmis = p.first;
        Long64_t nEvt = p.second;
        Double_t moyWater = nEvt > 0 ? sumWaterGivenEmitted[nEmis] / nEvt : 0.;
        Double_t taux = nEmis > 0 ? 100. * moyWater / nEmis : 0.;
        
        std::cout << "  │    " << std::setw(2) << nEmis << "   │ " << std::setw(8) << nEvt << "   │   "
                  << std::setw(6) << std::setprecision(3) << moyWater << "   │   "
                  << std::setw(6) << std::setprecision(2) << taux << "   │\n";
    }
    std::cout << "  └─────────┴────────────┴────────────┴────────────┘\n\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // ANALYSE 4: Matrice Raie × Multiplicité
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  4. MATRICE RAIE × MULTIPLICITÉ PAR ÉVÉNEMENT\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    // Construire la matrice
    const Int_t maxMult = maxGamma + 1;
    std::vector<std::vector<Long64_t>> matrixRaieMult(kNbGammaLines, std::vector<Long64_t>(maxMult, 0));
    
    for (auto& p : eventLineIDs) {
        Int_t evtID = p.first;
        Int_t mult = eventMultiplicity[evtID];
        
        if (mult < maxMult) {
            for (Int_t lineID : p.second) {
                if (lineID >= 0 && lineID < kNbGammaLines) {
                    matrixRaieMult[lineID][mult]++;
                }
            }
        }
    }
    
    // Affichage de la matrice
    std::cout << "  Nombre de gammas par raie en fonction de la multiplicité de l'événement:\n\n";
    std::cout << "  Raie\\Mult │";
    for (Int_t m = 0; m < maxMult; ++m) {
        std::cout << "  " << std::setw(3) << m << "  │";
    }
    std::cout << "\n  ──────────";
    for (Int_t m = 0; m < maxMult; ++m) std::cout << "────────";
    std::cout << "\n";
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        std::cout << "  " << std::setw(9) << kGammaLineNamesShort[l] << " │";
        for (Int_t m = 0; m < maxMult; ++m) {
            std::cout << " " << std::setw(5) << matrixRaieMult[l][m] << " │";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // ANALYSE 5: Gammas perdus dans l'air
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  5. GAMMAS AYANT INTERAGI DANS L'AIR AVANT D'ATTEINDRE L'EAU\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    Long64_t totalNotReachedWater = 0;
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        totalNotReachedWater += lineNotReachedWater[l];
    }
    
    std::cout << "  Nombre total de gammas émis:              " << totalGammas << "\n";
    std::cout << "  Gammas n'atteignant PAS l'eau:            " << totalNotReachedWater << "\n";
    std::cout << "  Fraction interagissant dans l'air:       " << std::setprecision(2) 
              << 100.*totalNotReachedWater/totalGammas << "%\n\n";
    
    // Détail par raie
    std::cout << "  Détail par raie:\n";
    std::cout << "  ┌──────┬─────────────┬─────────────┬───────────────┬────────────────┐\n";
    std::cout << "  │  ID  │ Raie        │ N_émis      │ N_non→eau     │ Taux perte (%) │\n";
    std::cout << "  ├──────┼─────────────┼─────────────┼───────────────┼────────────────┤\n";
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        Double_t tauxPerte = lineCounts[l] > 0 ? 100. * lineNotReachedWater[l] / lineCounts[l] : 0.;
        std::cout << "  │  " << std::setw(2) << l << "  │ " << std::setw(11) << std::left << kGammaLineNamesShort[l] << std::right
                  << " │ " << std::setw(9) << lineCounts[l] << "   │ "
                  << std::setw(11) << lineNotReachedWater[l] << "   │    "
                  << std::setw(8) << std::setprecision(2) << tauxPerte << "    │\n";
    }
    std::cout << "  └──────┴─────────────┴─────────────┴───────────────┴────────────────┘\n\n";
    
    // Distribution du nombre de gammas perdus par événement
    TH1D* hNLost = new TH1D("hNLost", "Gammas perdus dans l'air par evenement;N_{#gamma} perdus;Evenements",
                            10, -0.5, 9.5);
    hNLost->SetLineColor(kRed+1);
    hNLost->SetFillColor(kRed);
    hNLost->SetFillStyle(1001);
    
    std::map<Int_t, Long64_t> lostCounts;
    
    for (Long64_t i = 0; i < nEvents; ++i) {
        eventTree->GetEntry(i);
        Int_t nLost = evt_NGammaEmitted - evt_NGammaWater;
        hNLost->Fill(nLost);
        lostCounts[nLost]++;
    }
    
    std::cout << "  Distribution du nombre de gammas perdus dans l'air par événement:\n";
    std::cout << "  ┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "  │ N_perdu │  Nombre d'événements  │  Fraction (%)  │  Cumulé (%)     │\n";
    std::cout << "  ├─────────┼───────────────────────┼────────────────┼─────────────────┤\n";
    
    cumul = 0.;
    for (auto& p : lostCounts) {
        Double_t frac = 100. * p.second / nEvents;
        cumul += frac;
        std::cout << "  │    " << std::setw(2) << p.first << "   │       " 
                  << std::setw(8) << p.second << "        │     " 
                  << std::setw(6) << std::fixed << std::setprecision(2) << frac << "     │     "
                  << std::setw(6) << cumul << "      │\n";
    }
    std::cout << "  └─────────┴───────────────────────┴────────────────┴─────────────────┘\n\n";
    
    Double_t meanLost = hNLost->GetMean();
    std::cout << "  Moyenne de gammas perdus/événement: " << std::setprecision(4) << meanLost << "\n\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // GÉNÉRATION DES FIGURES
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  GÉNÉRATION DES FIGURES\n";
    std::cout << "════════════════════════════════════════════════════════════════════════\n\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // Figure 1: Distribution du nombre de gammas par événement
    // ───────────────────────────────────────────────────────────────────────
    
    TCanvas* c1 = new TCanvas("c1_emission", "Distribution gammas/evenement", 1400, 500);
    c1->Divide(2, 1);
    
    c1->cd(1);
    gPad->SetLogy();
    hNGammaEmitted->SetMinimum(0.5);
    hNGammaEmitted->Draw("BAR");
    
    TLine* lineMean = new TLine(meanGamma, 0, meanGamma, hNGammaEmitted->GetMaximum()*0.9);
    lineMean->SetLineColor(kRed);
    lineMean->SetLineStyle(2);
    lineMean->SetLineWidth(3);
    lineMean->Draw("SAME");
    
    TLine* lineTheo = new TLine(kMeanGammaTheory, 0, kMeanGammaTheory, hNGammaEmitted->GetMaximum()*0.9);
    lineTheo->SetLineColor(kMagenta+2);
    lineTheo->SetLineStyle(3);
    lineTheo->SetLineWidth(3);
    lineTheo->Draw("SAME");
    
    TLegend* leg1 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg1->AddEntry(hNGammaEmitted, "Simulation", "f");
    leg1->AddEntry(lineMean, Form("Moyenne simu = %.3f", meanGamma), "l");
    leg1->AddEntry(lineTheo, Form("Moyenne theo = %.3f", kMeanGammaTheory), "l");
    leg1->Draw();
    
    c1->cd(2);
    gPad->SetLogy();
    hNGammaEmitted->Draw("BAR");
    lineMean->Draw("SAME");
    lineTheo->Draw("SAME");
    
    c1->SaveAs("emission_1_gamma_per_event.png");
    std::cout << ">>> Sauvegardé: emission_1_gamma_per_event.png\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // Figure 2: Distribution par raie
    // ───────────────────────────────────────────────────────────────────────
    
    TCanvas* c2 = new TCanvas("c2_raies", "Distribution par raie", 1400, 500);
    c2->Divide(2, 1);
    
    // Histogramme des comptages absolus
    TH1D* hLineCounts = new TH1D("hLineCounts", "Gammas emis par raie;Raie;N_{#gamma}", 
                                  kNbGammaLines, -0.5, kNbGammaLines - 0.5);
    hLineCounts->SetLineColor(kOrange+2);
    hLineCounts->SetFillColor(kOrange);
    hLineCounts->SetFillStyle(1001);
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        hLineCounts->SetBinContent(l + 1, lineCounts[l]);
        hLineCounts->GetXaxis()->SetBinLabel(l + 1, kGammaLineNamesShort[l]);
    }
    
    c2->cd(1);
    gPad->SetLogy();
    hLineCounts->Draw("BAR");
    
    // Comparaison simulation/théorie
    c2->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    
    TH1D* hSimuPct = new TH1D("hSimuPct", "Comparaison Simulation vs Theorie;Raie;Intensite relative (%)",
                              kNbGammaLines, -0.5, kNbGammaLines - 0.5);
    TH1D* hTheoPct = new TH1D("hTheoPct", "", kNbGammaLines, -0.5, kNbGammaLines - 0.5);
    
    hSimuPct->SetLineColor(kBlue+1);
    hSimuPct->SetFillColor(kCyan);
    hSimuPct->SetFillStyle(1001);
    hTheoPct->SetLineColor(kRed+1);
    hTheoPct->SetFillColor(kMagenta-7);
    hTheoPct->SetFillStyle(1001);
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        Double_t simuPct = totalGammas > 0 ? 100. * lineCounts[l] / totalGammas : 0.;
        Double_t theoPct = 100. * kGammaLineIntensities[l] / sumTheoIntensity;
        
        hSimuPct->SetBinContent(l + 1, simuPct);
        hTheoPct->SetBinContent(l + 1, theoPct);
        hSimuPct->GetXaxis()->SetBinLabel(l + 1, kGammaLineNamesShort[l]);
    }
    
    hSimuPct->SetBarWidth(0.35);
    hSimuPct->SetBarOffset(0.15);
    hTheoPct->SetBarWidth(0.35);
    hTheoPct->SetBarOffset(0.52);
    
    hSimuPct->Draw("BAR");
    hTheoPct->Draw("BAR SAME");
    
    TLegend* leg2 = new TLegend(0.15, 0.75, 0.38, 0.88);
    leg2->AddEntry(hSimuPct, "Simulation", "f");
    leg2->AddEntry(hTheoPct, "Theorie LNHB", "f");
    leg2->Draw();
    
    c2->SaveAs("emission_2_gamma_per_line.png");
    std::cout << ">>> Sauvegardé: emission_2_gamma_per_line.png\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // Figure 3: Gammas atteignant l'eau
    // ───────────────────────────────────────────────────────────────────────
    
    TCanvas* c3 = new TCanvas("c3_water", "Gammas atteignant l'eau", 1400, 500);
    c3->Divide(2, 1);
    
    c3->cd(1);
    gPad->SetLogy();
    hNGammaWater->SetMinimum(0.5);
    hNGammaWater->Draw("BAR");
    
    TLine* lineMeanW = new TLine(meanWater, 0, meanWater, hNGammaWater->GetMaximum()*0.9);
    lineMeanW->SetLineColor(kRed);
    lineMeanW->SetLineStyle(2);
    lineMeanW->SetLineWidth(3);
    lineMeanW->Draw("SAME");
    
    TLegend* leg3 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg3->AddEntry(hNGammaWater, "Gammas #rightarrow eau", "f");
    leg3->AddEntry(lineMeanW, Form("Moyenne = %.3f", meanWater), "l");
    leg3->Draw();
    
    c3->cd(2);
    gPad->SetLogy();
    hNGammaEmitted->SetTitle("Comparaison: emis vs atteignant l'eau");
    hNGammaEmitted->Draw("BAR");
    hNGammaWater->Draw("BAR SAME");
    
    TLegend* leg3b = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg3b->AddEntry(hNGammaEmitted, "Emis", "f");
    leg3b->AddEntry(hNGammaWater, "#rightarrow Eau", "f");
    leg3b->Draw();
    
    c3->SaveAs("emission_3_gamma_reaching_water.png");
    std::cout << ">>> Sauvegardé: emission_3_gamma_reaching_water.png\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // Figure 4: Histogramme 2D Raie × Multiplicité
    // ───────────────────────────────────────────────────────────────────────
    
    TCanvas* c4 = new TCanvas("c4_2D", "Matrice Raie x Multiplicite", 1400, 600);
    c4->Divide(2, 1);
    
    TH2D* h2DRaieMult = new TH2D("h2DRaieMult", 
        "Matrice Raie #times Multiplicite;Multiplicite (N_{#gamma}/evt);Raie gamma",
        maxMult, -0.5, maxMult - 0.5, kNbGammaLines, -0.5, kNbGammaLines - 0.5);
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        h2DRaieMult->GetYaxis()->SetBinLabel(l + 1, kGammaLineNamesShort[l]);
        for (Int_t m = 0; m < maxMult; ++m) {
            h2DRaieMult->SetBinContent(m + 1, l + 1, matrixRaieMult[l][m]);
        }
    }
    
    c4->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    gStyle->SetPalette(kRainBow);
    h2DRaieMult->GetYaxis()->SetLabelSize(0.045);
    h2DRaieMult->GetYaxis()->SetLabelOffset(0.005);
    h2DRaieMult->GetZaxis()->SetLabelSize(0.03);
    h2DRaieMult->Draw("COLZ");
    
    // Ajouter les valeurs dans les cellules
    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(22);
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        for (Int_t m = 0; m < maxMult; ++m) {
            if (matrixRaieMult[l][m] > 0) {
                Double_t maxVal = h2DRaieMult->GetMaximum();
                if (matrixRaieMult[l][m] > maxVal/2) latex.SetTextColor(kWhite);
                else latex.SetTextColor(kBlack);
                latex.DrawLatex(m, l, Form("%lld", matrixRaieMult[l][m]));
            }
        }
    }
    
    c4->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    gStyle->SetPalette(kRainBow);
    TH2D* h2DRaieMultLog = (TH2D*)h2DRaieMult->Clone("h2DRaieMultLog");
    h2DRaieMultLog->SetTitle("Matrice Raie #times Multiplicite (log)");
    h2DRaieMultLog->GetYaxis()->SetLabelSize(0.045);
    h2DRaieMultLog->GetYaxis()->SetLabelOffset(0.005);
    h2DRaieMultLog->GetZaxis()->SetLabelSize(0.03);
    h2DRaieMultLog->Draw("COLZ");
    
    c4->SaveAs("emission_4_2D_line_vs_multiplicity.png");
    std::cout << ">>> Sauvegardé: emission_4_2D_line_vs_multiplicity.png\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // Figure 5: Gammas perdus dans l'air
    // ───────────────────────────────────────────────────────────────────────
    
    TCanvas* c5 = new TCanvas("c5_lost", "Gammas perdus dans l'air", 1400, 900);
    c5->Divide(2, 2);
    
    c5->cd(1);
    gPad->SetLogy();
    hNLost->SetMinimum(0.5);
    hNLost->Draw("BAR");
    TLine* lineMeanL = new TLine(meanLost, 0, meanLost, hNLost->GetMaximum()*0.9);
    lineMeanL->SetLineColor(kBlue+2);
    lineMeanL->SetLineStyle(2);
    lineMeanL->SetLineWidth(3);
    lineMeanL->Draw("SAME");
    
    TLegend* leg5 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg5->AddEntry(hNLost, "Gammas perdus/evt", "f");
    leg5->AddEntry(lineMeanL, Form("Moyenne = %.4f", meanLost), "l");
    leg5->Draw();
    
    // Taux de perte par raie
    c5->cd(2);
    gPad->SetLogy();
    TH1D* hLossRate = new TH1D("hLossRate", "Taux de perte par raie;Raie;Taux de perte (%)",
                               kNbGammaLines, -0.5, kNbGammaLines - 0.5);
    hLossRate->SetLineColor(kMagenta+1);
    hLossRate->SetFillColor(kMagenta);
    hLossRate->SetFillStyle(1001);
    hLossRate->SetMinimum(0.5);
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        Double_t rate = lineCounts[l] > 0 ? 100. * lineNotReachedWater[l] / lineCounts[l] : 0.;
        hLossRate->SetBinContent(l + 1, rate);
        hLossRate->GetXaxis()->SetBinLabel(l + 1, kGammaLineNamesShort[l]);
    }
    hLossRate->Draw("BAR");
    
    // Histogramme 2D émis vs perdus
    c5->cd(3);
    TH2D* h2DEmisLost = new TH2D("h2DEmisLost", 
        "N_{#gamma} emis vs N_{#gamma} perdus;N_{#gamma} emis;N_{#gamma} perdus",
        maxMult, -0.5, maxMult - 0.5, 5, -0.5, 4.5);
    
    for (Long64_t i = 0; i < nEvents; ++i) {
        eventTree->GetEntry(i);
        Int_t nLost = evt_NGammaEmitted - evt_NGammaWater;
        h2DEmisLost->Fill(evt_NGammaEmitted, nLost);
    }
    
    gPad->SetLogz();
    gStyle->SetPalette(kRainBow);
    h2DEmisLost->Draw("COLZ");
    
    // Moyenne conditionnelle
    c5->cd(4);
    TGraphErrors* gMeanLost = new TGraphErrors();
    gMeanLost->SetTitle("Moyenne gammas perdus vs emis;N_{#gamma} emis;Moyenne perdus");
    gMeanLost->SetMarkerStyle(20);
    gMeanLost->SetMarkerSize(1.8);
    gMeanLost->SetMarkerColor(kBlue+2);
    gMeanLost->SetLineColor(kBlue+2);
    gMeanLost->SetLineWidth(3);
    
    Int_t iPoint = 0;
    for (auto& p : multCounts) {
        Int_t nEmis = p.first;
        if (nEmis > 0 && p.second > 0) {
            // Calculer moyenne et erreur des perdus pour cette multiplicité
            Double_t sumL = 0., sumL2 = 0.;
            Long64_t countL = 0;
            
            for (Long64_t i = 0; i < nEvents; ++i) {
                eventTree->GetEntry(i);
                if (evt_NGammaEmitted == nEmis) {
                    Int_t nLost = evt_NGammaEmitted - evt_NGammaWater;
                    sumL += nLost;
                    sumL2 += nLost * nLost;
                    countL++;
                }
            }
            
            if (countL > 0) {
                Double_t meanL = sumL / countL;
                Double_t varL = sumL2 / countL - meanL * meanL;
                Double_t errL = varL > 0 ? TMath::Sqrt(varL / countL) : 0.;
                
                gMeanLost->SetPoint(iPoint, nEmis, meanL);
                gMeanLost->SetPointError(iPoint, 0, errL);
                iPoint++;
            }
        }
    }
    
    gMeanLost->Draw("APE");
    
    c5->SaveAs("emission_5_gamma_lost_in_air.png");
    std::cout << ">>> Sauvegardé: emission_5_gamma_lost_in_air.png\n";
    
    // ───────────────────────────────────────────────────────────────────────
    // Figure 6: Résumé synthétique
    // ───────────────────────────────────────────────────────────────────────
    
    TCanvas* c6 = new TCanvas("c6_summary", "Resume emission Am-241", 1400, 900);
    c6->Divide(2, 2);
    
    // Pie chart des raies
    c6->cd(1);
    
    const Int_t nPie = kNbGammaLines;
    Double_t vals[nPie];
    const char* labels[nPie];
    Int_t colors[nPie] = {kCyan, kCyan+1, kCyan+2, kCyan+3, kGreen, kGreen+2, 
                          kYellow, kYellow+2, kRed, kOrange, kOrange+2, kMagenta};
    
    for (Int_t l = 0; l < nPie; ++l) {
        vals[l] = lineCounts[l];
        labels[l] = Form("%.1f keV", kGammaLineEnergies[l]);
    }
    
    TPie* pie = new TPie("pie", "Repartition des raies emises", nPie, vals, colors);
    pie->SetLabelsOffset(0.02);
    pie->SetLabelFormat("%txt (%perc)");
    for (Int_t l = 0; l < nPie; ++l) {
        pie->SetEntryLabel(l, labels[l]);
        if (vals[l] / totalGammas < 0.02) pie->SetEntryLabel(l, ""); // Masquer les petites
    }
    //pie->Draw("3d");
    
    // Tableau récapitulatif
    c6->cd(2);
    TPaveText* pt = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    
    pt->AddText(Form("Nombre d'evenements:    %lld", nEvents));
    pt->AddText(Form("Gammas emis (total):    %lld", totalGammas));
    pt->AddText(Form("Moyenne gamma/evt:       %.4f", meanGamma));
    pt->AddText(Form("Theorie gamma/evt:       %.4f", kMeanGammaTheory));
    pt->AddText(Form("Ecart relatif:          %.2f%%", 100.*(meanGamma-kMeanGammaTheory)/kMeanGammaTheory));
    pt->AddText("");
    pt->AddText(Form("Evt sans gamma:         %lld (%.1f%%)", multCounts[0], 100.*multCounts[0]/nEvents));
    pt->AddText(Form("Gamma atteignant eau:   %lld", (Long64_t)sumWater));
    pt->AddText(Form("Taux passage -> eau:    %.1f%%", 100.*sumWater/totalGammas));
    pt->AddText(Form("Gamma perdus (air):     %lld (%.1f%%)", totalNotReachedWater, 100.*totalNotReachedWater/totalGammas));
    
    //pt->Draw();
    
    // Distribution avec fit Poisson
    c6->cd(3);
    gPad->SetLogy(0);
    TH1D* hNGammaNorm = (TH1D*)hNGammaEmitted->Clone("hNGammaNorm");
    hNGammaNorm->SetTitle("Distribution vs Poisson;N_{#gamma};Densite");
    hNGammaNorm->Scale(1./nEvents);
    hNGammaNorm->SetMinimum(1e-6); // Pour éviter les problèmes avec log(0)
    hNGammaNorm->Draw("BAR");
    
    // Poisson
    TF1* fPoisson = new TF1("fPoisson", "[0]*TMath::Poisson(x,[1])", -0.5, maxGamma + 0.5);
    fPoisson->SetParameters(1., meanGamma);
    fPoisson->SetLineColor(kRed+1);
    fPoisson->SetLineWidth(3);
    fPoisson->Draw("SAME");
    
    TLegend* leg6 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg6->AddEntry(hNGammaNorm, "Simulation", "f");
    leg6->AddEntry(fPoisson, Form("Poisson(#lambda=%.3f)", meanGamma), "l");
    leg6->Draw();
    
    // Énergies des raies
    c6->cd(4);
    TH1D* hEnergies = new TH1D("hEnergies", "Energies des raies Am-241;Raie;Energie (keV)",
                               kNbGammaLines, -0.5, kNbGammaLines - 0.5);
    hEnergies->SetLineColor(kOrange+1);
    hEnergies->SetFillColor(kYellow);
    hEnergies->SetFillStyle(1001);
    
    for (Int_t l = 0; l < kNbGammaLines; ++l) {
        hEnergies->SetBinContent(l + 1, kGammaLineEnergies[l]);
        hEnergies->GetXaxis()->SetBinLabel(l + 1, kGammaLineNamesShort[l]);
    }
    hEnergies->Draw("BAR");
    
    TLine* line59 = new TLine(-0.5, 59.5409, kNbGammaLines - 0.5, 59.5409);
    line59->SetLineColor(kRed);
    line59->SetLineStyle(2);
    line59->SetLineWidth(3);
    line59->Draw("SAME");
    
    TLegend* leg6b = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg6b->AddEntry(line59, "Raie principale 59.5 keV", "l");
    leg6b->Draw();
    
    c6->SaveAs("emission_6_summary.png");
    std::cout << ">>> Sauvegardé: emission_6_summary.png\n";
    
    // ═══════════════════════════════════════════════════════════════════════
    // RÉSUMÉ FINAL
    // ═══════════════════════════════════════════════════════════════════════
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║     ANALYSE TERMINÉE AVEC SUCCÈS                                     ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Fichiers générés:                                                   ║\n";
    std::cout << "║    - emission_1_gamma_per_event.png     : Distrib. γ/événement       ║\n";
    std::cout << "║    - emission_2_gamma_per_line.png      : Distrib. γ par raie        ║\n";
    std::cout << "║    - emission_3_gamma_reaching_water.png: γ atteignant l'eau         ║\n";
    std::cout << "║    - emission_4_2D_line_vs_multiplicity.png: Matrice 2D raie×mult    ║\n";
    std::cout << "║    - emission_5_gamma_lost_in_air.png   : γ perdus dans l'air        ║\n";
    std::cout << "║    - emission_6_summary.png             : Résumé synthétique         ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════╝\n\n";
}
