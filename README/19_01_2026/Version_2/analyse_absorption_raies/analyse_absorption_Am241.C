/*******************************************************************************
 * analyse_absorption_Am241.C
 * 
 * Script ROOT pour l'analyse des résultats de simulation Geant4 Am-241
 * VERSION CORRIGÉE avec:
 *   - Visualisation recommandée (barres par raie au lieu d'histogrammes)
 *   - Explication du problème de pollution Compton
 *   - Comparaison des deux méthodes de visualisation
 * 
 * Usage: root -l 'analyse_absorption_Am241_v2.C("output.root")'
 * 
 ******************************************************************************/

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TColor.h"
#include "TROOT.h"
#include "TMath.h"
#include "TBox.h"
#include "TArrow.h"
#include "TGaxis.h"
#include <iostream>
#include <iomanip>

// ═══════════════════════════════════════════════════════════════════════════
// DÉFINITION DES RAIES Am-241
// ═══════════════════════════════════════════════════════════════════════════
const int kNbRaies = 12;
const double kRaieEnergies[kNbRaies] = {
    11.89, 13.9, 17.0, 20.8, 26.3446, 33.1963, 
    43.420, 55.56, 59.5409, 98.97, 102.98, 125.30
};
const char* kRaieNames[kNbRaies] = {
    "X_{Ll}", "X_{L#alpha}", "X_{L#beta}", "X_{L#gamma}",
    "#gamma 26", "#gamma 33", "#gamma 43", "#gamma 56",
    "#gamma 59.5", "#gamma 99", "#gamma 103", "#gamma 125"
};
const char* kRaieNamesShort[kNbRaies] = {
    "X_Ll", "X_La", "X_Lb", "X_Lg",
    "g26", "g33", "g43", "g56",
    "g59", "g99", "g103", "g125"
};

// Coefficients d'atténuation massique de l'eau (cm²/g) - NIST
const double kMuRho[kNbRaies] = {
    4.0, 2.5, 1.4, 0.80, 0.45, 0.28, 0.22, 0.20, 0.19, 0.17, 0.17, 0.16
};

// ═══════════════════════════════════════════════════════════════════════════
// FONCTION PRINCIPALE
// ═══════════════════════════════════════════════════════════════════════════
void analyse_absorption_Am241(const char* filename = "output.root") {
    
    // ─────────────────────────────────────────────────────────────────────────
    // Configuration du style
    // ─────────────────────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetLabelSize(0.035, "XY");
    gStyle->SetTitleSize(0.04, "XY");
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetHistLineWidth(2);
    gStyle->SetPalette(kRainBow);
    
    // Couleurs
    Int_t colEmis = TColor::GetColor("#3498db");
    Int_t colEntres = TColor::GetColor("#2ecc71");
    Int_t colAbsorbes = TColor::GetColor("#e74c3c");
    Int_t colTaux = TColor::GetColor("#9b59b6");
    Int_t colTheo = TColor::GetColor("#f39c12");
    Int_t colCompton = TColor::GetColor("#e67e22");
    
    // ─────────────────────────────────────────────────────────────────────────
    // Ouverture du fichier
    // ─────────────────────────────────────────────────────────────────────────
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "ERREUR: Impossible d'ouvrir " << filename << std::endl;
        return;
    }
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "  ANALYSE DES TAUX D'ABSORPTION - Am-241 (VERSION CORRIGEE)" << std::endl;
    std::cout << "  Fichier: " << filename << std::endl;
    std::cout << std::string(80, '=') << "\n" << std::endl;
    
    // ─────────────────────────────────────────────────────────────────────────
    // Récupération des données
    // ─────────────────────────────────────────────────────────────────────────
    TTree* tGammaLines = (TTree*)file->Get("gamma_lines");
    TTree* tEventData = (TTree*)file->Get("EventData");
    TH1D* hGammaEmitted = (TH1D*)file->Get("hGammaEmitted");
    TH1D* hGammaEnteringWater = (TH1D*)file->Get("hGammaEnteringWater");
    
    if (!tGammaLines || !tEventData) {
        std::cerr << "ERREUR: TTrees manquants" << std::endl;
        return;
    }
    
    Long64_t nEvents = tEventData->GetEntries();
    std::cout << "Nombre d'evenements: " << nEvents << std::endl;
    
    // Lecture des données gamma_lines
    Int_t lineIndex, emitted, enteredWater, absorbedWater;
    Double_t energy_keV, waterAbsRate;
    
    tGammaLines->SetBranchAddress("lineIndex", &lineIndex);
    tGammaLines->SetBranchAddress("energy_keV", &energy_keV);
    tGammaLines->SetBranchAddress("emitted", &emitted);
    tGammaLines->SetBranchAddress("enteredWater", &enteredWater);
    tGammaLines->SetBranchAddress("absorbedWater", &absorbedWater);
    tGammaLines->SetBranchAddress("waterAbsRate", &waterAbsRate);
    
    // Tableaux pour stocker les données
    Double_t aEmis[kNbRaies], aEntres[kNbRaies], aAbsorbes[kNbRaies], aTaux[kNbRaies], aTauxTheo[kNbRaies];
    
    // ─────────────────────────────────────────────────────────────────────────
    // Création des histogrammes par raie
    // ─────────────────────────────────────────────────────────────────────────
    TH1D* hEmisParRaie = new TH1D("hEmisParRaie", 
        "Gammas emis par raie;Raie;Nombre", kNbRaies, -0.5, kNbRaies-0.5);
    TH1D* hEntresParRaie = new TH1D("hEntresParRaie", 
        "Gammas entres par raie;Raie;Nombre", kNbRaies, -0.5, kNbRaies-0.5);
    TH1D* hAbsorbesParRaie = new TH1D("hAbsorbesParRaie", 
        "Gammas absorbes par raie;Raie;Nombre", kNbRaies, -0.5, kNbRaies-0.5);
    TH1D* hTauxMesure = new TH1D("hTauxMesure", 
        "Taux mesure;Raie;Taux (%)", kNbRaies, -0.5, kNbRaies-0.5);
    TH1D* hTauxTheo = new TH1D("hTauxTheo", 
        "Taux theorique;Raie;Taux (%)", kNbRaies, -0.5, kNbRaies-0.5);
    
    // Labels des axes
    for (int i = 0; i < kNbRaies; i++) {
        hEmisParRaie->GetXaxis()->SetBinLabel(i+1, kRaieNames[i]);
        hEntresParRaie->GetXaxis()->SetBinLabel(i+1, kRaieNames[i]);
        hAbsorbesParRaie->GetXaxis()->SetBinLabel(i+1, kRaieNames[i]);
        hTauxMesure->GetXaxis()->SetBinLabel(i+1, kRaieNames[i]);
        hTauxTheo->GetXaxis()->SetBinLabel(i+1, kRaieNames[i]);
    }
    
    // ─────────────────────────────────────────────────────────────────────────
    // Remplissage des histogrammes et calcul des taux théoriques
    // ─────────────────────────────────────────────────────────────────────────
    std::cout << "\n" << std::string(95, '-') << std::endl;
    std::cout << std::setw(6) << "Raie" << std::setw(10) << "E (keV)" 
              << std::setw(10) << "Emis" << std::setw(10) << "Entres" 
              << std::setw(10) << "Absorbes" << std::setw(12) << "Mesure(%)" 
              << std::setw(12) << "Theo 3mm(%)" << std::setw(10) << "Ratio" 
              << std::setw(12) << "Statut" << std::endl;
    std::cout << std::string(95, '-') << std::endl;
    
    double totalEntres = 0, totalAbsorbes = 0;
    
    for (Long64_t i = 0; i < tGammaLines->GetEntries(); i++) {
        tGammaLines->GetEntry(i);
        
        // Stockage
        aEmis[lineIndex] = emitted;
        aEntres[lineIndex] = enteredWater;
        aAbsorbes[lineIndex] = absorbedWater;
        aTaux[lineIndex] = waterAbsRate;
        
        // Calcul taux théorique (3 mm eau)
        double epaisseur_cm = 0.3;
        double tauxTheo = 100.0 * (1.0 - TMath::Exp(-kMuRho[lineIndex] * epaisseur_cm));
        aTauxTheo[lineIndex] = tauxTheo;
        
        // Remplissage histogrammes
        hEmisParRaie->SetBinContent(lineIndex+1, emitted);
        hEntresParRaie->SetBinContent(lineIndex+1, enteredWater);
        hAbsorbesParRaie->SetBinContent(lineIndex+1, absorbedWater > 0 ? absorbedWater : 0.1);
        hTauxMesure->SetBinContent(lineIndex+1, waterAbsRate);
        hTauxTheo->SetBinContent(lineIndex+1, tauxTheo);
        
        // Ratio et statut
        double ratio = (tauxTheo > 0) ? waterAbsRate / tauxTheo : 0;
        totalEntres += enteredWater;
        totalAbsorbes += absorbedWater;
        
        const char* statut = "";
        if (absorbedWater == 0) statut = "Stat. insuf.";
        else if (absorbedWater < 10) statut = "Faible stat.";
        else statut = "OK";
        
        std::cout << std::setw(6) << lineIndex 
                  << std::setw(10) << std::fixed << std::setprecision(1) << energy_keV
                  << std::setw(10) << emitted << std::setw(10) << enteredWater 
                  << std::setw(10) << absorbedWater 
                  << std::setw(12) << std::setprecision(2) << waterAbsRate 
                  << std::setw(12) << tauxTheo
                  << std::setw(10) << std::setprecision(2) << ratio 
                  << std::setw(12) << statut << std::endl;
    }
    
    std::cout << std::string(95, '-') << std::endl;
    std::cout << std::setw(6) << "TOTAL" << std::setw(10) << "-"
              << std::setw(10) << "-" << std::setw(10) << (int)totalEntres 
              << std::setw(10) << (int)totalAbsorbes 
              << std::setw(12) << std::setprecision(2) << 100.0*totalAbsorbes/totalEntres 
              << std::endl;
    
    // ─────────────────────────────────────────────────────────────────────────
    // Calcul du maximum pour les échelles Y (adaptatif)
    // ─────────────────────────────────────────────────────────────────────────
    double maxVal = 0;
    for (int i = 0; i < kNbRaies; i++) {
        if (aEmis[i] > maxVal) maxVal = aEmis[i];
        if (aEntres[i] > maxVal) maxVal = aEntres[i];
    }
    double yMax = maxVal * 5;      // Pour les barres avec annotations
    double yMaxHisto = maxVal * 2; // Pour les histogrammes sans annotations
    
    std::cout << "\n  Max gammas par raie: " << (int)maxVal << " -> yMax = " << (int)yMax << std::endl;
    
    // ═════════════════════════════════════════════════════════════════════════
    // CANVAS 1: PROBLÈME DE VISUALISATION - Histogramme vs Barres par raie
    // ═════════════════════════════════════════════════════════════════════════
    
    TCanvas* c1 = new TCanvas("c1", "Probleme de visualisation", 1400, 900);
    c1->Divide(2, 2);
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 1: Spectre émis (histogramme en énergie) - zoom 40-75 keV
    // ─────────────────────────────────────────────────────────────────────────
    c1->cd(1);
    gPad->SetLogy();
    gPad->SetGridy();
    
    TH1D* hEmisZoom = (TH1D*)hGammaEmitted->Clone("hEmisZoom");
    hEmisZoom->SetTitle("Spectre EMIS (histogramme en energie);Energie (keV);Nombre de gammas");
    hEmisZoom->GetXaxis()->SetRangeUser(40, 75);
    hEmisZoom->SetFillColor(colEmis);
    hEmisZoom->SetFillStyle(1001);
    hEmisZoom->SetMinimum(1);
    hEmisZoom->SetMaximum(yMaxHisto);
    hEmisZoom->Draw("BAR");
    
    // Lignes verticales pour les raies
    TLine* l1 = new TLine(43.4, 1, 43.4, yMaxHisto/10); l1->SetLineColor(kRed); l1->SetLineStyle(2); l1->Draw();
    TLine* l2 = new TLine(55.6, 1, 55.6, yMaxHisto/10); l2->SetLineColor(kRed); l2->SetLineStyle(2); l2->Draw();
    TLine* l3 = new TLine(59.5, 1, 59.5, yMaxHisto); l3->SetLineColor(kRed); l3->SetLineStyle(2); l3->Draw();
    
    TLatex latex;
    latex.SetTextSize(0.045);
    latex.SetTextColor(kRed);
    latex.DrawLatex(43.4, 2e4, "43.4");
    latex.DrawLatex(55.6, 2e4, "55.6");
    latex.DrawLatex(59.5, 5e4, "59.5");
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 2: Spectre entré dans l'eau - AVEC POLLUTION COMPTON
    // ─────────────────────────────────────────────────────────────────────────
    c1->cd(2);
    gPad->SetLogy();
    gPad->SetGridy();
    
    TH1D* hWaterZoom = (TH1D*)hGammaEnteringWater->Clone("hWaterZoom");
    hWaterZoom->SetTitle("Spectre ENTRE dans l'eau - #color[2]{POLLUE par Compton!};Energie (keV);Nombre de gammas");
    hWaterZoom->GetXaxis()->SetRangeUser(40, 75);
    hWaterZoom->SetFillColor(colAbsorbes);
    hWaterZoom->SetFillStyle(1001);
    hWaterZoom->SetMinimum(1);
    hWaterZoom->SetMaximum(yMaxHisto);
    hWaterZoom->Draw("BAR");
    
    // Zone Compton
    TBox* boxCompton = new TBox(48, 1, 58, yMaxHisto/10);
    boxCompton->SetFillColor(colCompton);
    boxCompton->SetFillStyle(3004);
    boxCompton->Draw();
    
    // Lignes verticales
    l1->Draw(); l2->Draw(); l3->Draw();
    
    // Légende
    TLegend* leg1 = new TLegend(0.12, 0.75, 0.52, 0.88);
    leg1->SetBorderSize(1);
    leg1->SetFillColor(kWhite);
    leg1->AddEntry(boxCompton, "Fond Compton (gammas 59.5 keV diffuses)", "f");
    leg1->Draw();
    
    latex.SetTextColor(colCompton);
    latex.DrawLatex(50, yMaxHisto/95, "Compton!");
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 3: SOLUTION - Barres par raie (entrés)
    // ─────────────────────────────────────────────────────────────────────────
    c1->cd(3);
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.20);
    
    hEntresParRaie->SetTitle("#checkmark CORRECT: Entres par RAIE (gamma_lines);Raie;Nombre de gammas");
    hEntresParRaie->SetFillColor(colEntres);
    hEntresParRaie->SetFillStyle(1001);
    hEntresParRaie->SetMinimum(0.5);
    hEntresParRaie->SetMaximum(yMax);
    hEntresParRaie->GetXaxis()->SetLabelSize(0.06);
    hEntresParRaie->GetXaxis()->LabelsOption("v");
    hEntresParRaie->Draw("BAR");
    
    // Valeurs sur les barres
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.028);
    for (int i = 0; i < kNbRaies; i++) {
        if (aEntres[i] > 0) {
            latex.DrawLatex(i - 0.3, aEntres[i]*1.5, Form("%.0f", aEntres[i]));
        }
    }
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 4: SOLUTION - Barres par raie (absorbés)
    // ─────────────────────────────────────────────────────────────────────────
    c1->cd(4);
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetBottomMargin(0.20);
    
    hAbsorbesParRaie->SetTitle("#checkmark CORRECT: Absorbes par RAIE (gamma_lines);Raie;Nombre de gammas");
    hAbsorbesParRaie->SetFillColor(colAbsorbes);
    hAbsorbesParRaie->SetFillStyle(1001);
    hAbsorbesParRaie->SetMinimum(0.05);
    hAbsorbesParRaie->SetMaximum(yMax);
    hAbsorbesParRaie->GetXaxis()->SetLabelSize(0.06);
    hAbsorbesParRaie->GetXaxis()->LabelsOption("v");
    hAbsorbesParRaie->Draw("BAR");
    
    // Valeurs sur les barres
    latex.SetTextSize(0.028);
    for (int i = 0; i < kNbRaies; i++) {
        latex.DrawLatex(i - 0.3, TMath::Max(aAbsorbes[i], 0.5)*2, Form("%.0f", aAbsorbes[i]));
    }
    
    c1->SaveAs("absorption_probleme_visualisation.png");
    
    // ═════════════════════════════════════════════════════════════════════════
    // CANVAS 2: EXPLICATION DU PROBLÈME COMPTON
    // ═════════════════════════════════════════════════════════════════════════
    
    TCanvas* c2 = new TCanvas("c2", "Explication Compton", 1400, 900);
    c2->Divide(3, 2);
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 1: Comparaison émis vs entrés (zoom)
    // ─────────────────────────────────────────────────────────────────────────
    c2->cd(1);
    gPad->SetLogy();
    gPad->SetGridy();
    
    TH1D* hEmisZoom2 = (TH1D*)hGammaEmitted->Clone("hEmisZoom2");
    TH1D* hWaterZoom2 = (TH1D*)hGammaEnteringWater->Clone("hWaterZoom2");
    
    hEmisZoom2->SetTitle("Comparaison Emis vs Entres (45-68 keV);Energie (keV);Nombre de gammas");
    hEmisZoom2->GetXaxis()->SetRangeUser(45, 68);
    hEmisZoom2->SetFillColor(colEmis);
    hEmisZoom2->SetFillStyle(3001);
    hEmisZoom2->SetMinimum(1);
    hEmisZoom2->SetMaximum(yMaxHisto * 20);  // Facteur 20 pour réduire la hauteur visuelle à ~50%
    hEmisZoom2->Draw("BAR");
    
    hWaterZoom2->GetXaxis()->SetRangeUser(45, 68);
    hWaterZoom2->SetFillColor(colAbsorbes);
    hWaterZoom2->SetFillStyle(3004);
    hWaterZoom2->Draw("BAR SAME");
    
    TLegend* leg2 = new TLegend(0.15, 0.55, 0.45, 0.68);
    leg2->AddEntry(hEmisZoom2, "Emis", "f");
    leg2->AddEntry(hWaterZoom2, "Entres eau", "f");
    leg2->Draw();
    
    // Annotations
    TArrow* arr1 = new TArrow(53, yMaxHisto*5, 53, yMaxHisto/2, 0.02, "|>");
    arr1->SetLineColor(colCompton);
    arr1->SetLineWidth(2);
    arr1->Draw();
    latex.SetTextColor(colCompton);
    latex.SetTextSize(0.045);
    latex.DrawLatex(48, yMaxHisto*8, "Exces Compton!");
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 2: Différence (Entrés - Émis)
    // ─────────────────────────────────────────────────────────────────────────
    c2->cd(2);
    gPad->SetGridy();
    
    // Créer histogramme de différence
    TH1D* hDiff = (TH1D*)hGammaEnteringWater->Clone("hDiff");
    hDiff->Add(hGammaEmitted, -1);
    hDiff->SetTitle("Difference (Entres - Emis);Energie (keV);Difference");
    hDiff->GetXaxis()->SetRangeUser(45, 68);
    
    // Calculer min/max et appliquer facteur pour réduire la hauteur à ~50%
    double diffMin = hDiff->GetMinimum();
    double diffMax = hDiff->GetMaximum();
    double diffRange = diffMax - diffMin;
    hDiff->SetMinimum(diffMin * 2.5);  // Étendre vers le bas
    hDiff->SetMaximum(diffMax * 5);    // Étendre vers le haut
    
    // Colorer selon le signe
    hDiff->SetFillColor(colCompton);
    hDiff->SetLineColor(kBlack);
    hDiff->Draw("BAR");
    
    // Ligne à zéro
    TLine* lZero = new TLine(45, 0, 68, 0);
    lZero->SetLineColor(kBlack);
    lZero->SetLineWidth(2);
    lZero->Draw();
    
    // Annotations
    latex.SetTextColor(colCompton);
    latex.DrawLatex(50, diffMax * 1.2, "Fond Compton");
    latex.DrawLatex(50, diffMin * 0.3, "(exces)");
    latex.SetTextColor(colEmis);
    latex.DrawLatex(58, diffMin * 1.5, "Perte");
    latex.DrawLatex(58, diffMin * 1.8, "(deficit)");
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 3: Courbe Compton théorique
    // ─────────────────────────────────────────────────────────────────────────
    c2->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    
    // Calcul énergie Compton vs angle
    const int nPoints = 100;
    Double_t theta[nPoints], Ecompton[nPoints];
    Double_t E0 = 59.5;  // keV
    Double_t me = 511.0; // keV
    
    for (int i = 0; i < nPoints; i++) {
        theta[i] = i * 180.0 / (nPoints - 1);
        Double_t cosTheta = TMath::Cos(theta[i] * TMath::Pi() / 180.0);
        Ecompton[i] = E0 / (1.0 + (E0/me) * (1.0 - cosTheta));
    }
    
    TGraph* gCompton = new TGraph(nPoints, theta, Ecompton);
    gCompton->SetTitle("Energie Compton vs angle de diffusion;Angle #theta (deg);Energie apres Compton (keV)");
    gCompton->SetLineColor(colEmis);
    gCompton->SetLineWidth(3);
    gCompton->GetXaxis()->SetRangeUser(0, 180);
    gCompton->GetYaxis()->SetRangeUser(30, 74);  // Plage élargie pour réduire la hauteur visuelle
    gCompton->Draw("AL");
    
    // Lignes horizontales pour les raies
    TLine* lRaie59 = new TLine(0, 59.5, 180, 59.5);
    lRaie59->SetLineColor(colEntres); lRaie59->SetLineStyle(2); lRaie59->SetLineWidth(2); lRaie59->Draw();
    TLine* lRaie56 = new TLine(0, 55.6, 180, 55.6);
    lRaie56->SetLineColor(colCompton); lRaie56->SetLineStyle(2); lRaie56->SetLineWidth(2); lRaie56->Draw();
    TLine* lRaie43 = new TLine(0, 43.4, 180, 43.4);
    lRaie43->SetLineColor(colAbsorbes); lRaie43->SetLineStyle(2); lRaie43->SetLineWidth(2); lRaie43->Draw();
    
    // Zone de pollution
    TBox* boxPollution = new TBox(0, 48, 180, 58);
    boxPollution->SetFillColor(colCompton);
    boxPollution->SetFillStyle(3004);
    boxPollution->Draw();
    
    TLegend* leg3 = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg3->AddEntry(gCompton, "E_{Compton}(#theta)", "l");
    leg3->AddEntry(lRaie59, "Raie 59.5 keV", "l");
    leg3->AddEntry(lRaie56, "Raie 55.6 keV", "l");
    leg3->AddEntry(boxPollution, "Zone pollution", "f");
    leg3->Draw();
    
    // Annotations
    latex.SetTextSize(0.04);
    latex.SetTextColor(kBlack);
    latex.DrawLatex(60, 56.5, "#theta=60#circ #rightarrow 56 keV");
    latex.DrawLatex(90, 53.5, "#theta=90#circ #rightarrow 53 keV");
    
    c2->SaveAs("absorption_explication_compton.png");
    
    // ═════════════════════════════════════════════════════════════════════════
    // CANVAS 3: VISUALISATION RECOMMANDÉE - Taux d'absorption par raie
    // ═════════════════════════════════════════════════════════════════════════
    
    TCanvas* c3 = new TCanvas("c3", "Visualisation recommandee", 1400, 900);
    c3->Divide(1, 2);
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 1: Taux d'absorption par raie
    // ─────────────────────────────────────────────────────────────────────────
    c3->cd(1);
    gPad->SetBottomMargin(0.18);
    gPad->SetGridy();
    
    // Créer histogramme avec couleurs selon le taux
    TH1D* hTauxCouleur = (TH1D*)hTauxMesure->Clone("hTauxCouleur");
    hTauxCouleur->SetTitle("TAUX D'ABSORPTION par raie Am-241;Raie;Taux d'absorption (%)");
    hTauxCouleur->SetFillColor(colTaux);
    hTauxCouleur->SetBarWidth(0.8);
    hTauxCouleur->SetBarOffset(0.1);
    hTauxCouleur->SetMinimum(0);
    hTauxCouleur->SetMaximum(70);
    hTauxCouleur->GetXaxis()->SetLabelSize(0.06);
    hTauxCouleur->GetXaxis()->LabelsOption("v");
    hTauxCouleur->Draw("BAR");
    
    // Superposer taux théorique
    hTauxTheo->SetFillColor(0);
    hTauxTheo->SetFillStyle(0);
    hTauxTheo->SetLineColor(colTheo);
    hTauxTheo->SetLineWidth(3);
    hTauxTheo->SetLineStyle(2);
    hTauxTheo->Draw("HIST SAME");
    
    // Légende
    TLegend* leg4 = new TLegend(0.65, 0.75, 0.92, 0.88);
    leg4->AddEntry(hTauxCouleur, "Mesure", "f");
    leg4->AddEntry(hTauxTheo, "Theorique (3mm eau)", "l");
    leg4->Draw();
    
    // Valeurs sur les barres
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.025);
    for (int i = 0; i < kNbRaies; i++) {
        if (aTaux[i] > 0) {
            latex.DrawLatex(i - 0.4, aTaux[i] + 2, Form("%.1f%%", aTaux[i]));
            latex.SetTextSize(0.028);
            latex.DrawLatex(i - 0.4, aTaux[i] + 5, Form("(%d)", (int)aAbsorbes[i]));
            latex.SetTextSize(0.025);
        } else {
            latex.SetTextColor(kGray+1);
            latex.DrawLatex(i - 0.3, 2, "0");
            latex.SetTextColor(kBlack);
        }
    }
    
    // ─────────────────────────────────────────────────────────────────────────
    // Pad 2: Comparaison Émis / Entrés / Absorbés (avec annotations)
    // ─────────────────────────────────────────────────────────────────────────
    c3->cd(2);
    gPad->SetLogy();
    gPad->SetBottomMargin(0.18);
    gPad->SetGridy();
    
    // Décaler les histogrammes pour les grouper
    hEmisParRaie->SetBarWidth(0.25);
    hEmisParRaie->SetBarOffset(0.1);
    hEmisParRaie->SetFillColor(colEmis);
    hEmisParRaie->SetTitle("Comparaison: Emis #rightarrow Entres #rightarrow Absorbes;Raie;Nombre de gammas (echelle log)");
    hEmisParRaie->SetMinimum(0.3);
    hEmisParRaie->SetMaximum(yMax);
    hEmisParRaie->GetXaxis()->SetLabelSize(0.06);
    hEmisParRaie->GetXaxis()->LabelsOption("v");
    hEmisParRaie->Draw("BAR");
    
    hEntresParRaie->SetBarWidth(0.25);
    hEntresParRaie->SetBarOffset(0.35);
    hEntresParRaie->SetFillColor(colEntres);
    hEntresParRaie->Draw("BAR SAME");
    
    hAbsorbesParRaie->SetBarWidth(0.25);
    hAbsorbesParRaie->SetBarOffset(0.6);
    hAbsorbesParRaie->SetFillColor(colAbsorbes);
    hAbsorbesParRaie->Draw("BAR SAME");
    
    // *** ANNOTATIONS: Valeurs numériques sur chaque barre ***
    latex.SetTextSize(0.020);
    for (int i = 0; i < kNbRaies; i++) {
        // Position X des barres (avec offset)
        double xEmis = i - 0.5 + 0.1 + 0.125;      // centre de la barre Emis
        double xEntres = i - 0.5 + 0.35 + 0.125;   // centre de la barre Entres
        double xAbsorbes = i - 0.5 + 0.6 + 0.125;  // centre de la barre Absorbes
        
        // Émis (bleu)
        latex.SetTextColor(colEmis);
        double yEmis = aEmis[i] * 1.5;
        latex.DrawLatex(xEmis - 0.1, yEmis, Form("%.0f", aEmis[i]));
        
        // Entrés (vert)
        latex.SetTextColor(colEntres);
        double yEntres = aEntres[i] * 1.5;
        latex.DrawLatex(xEntres - 0.1, yEntres, Form("%.0f", aEntres[i]));
        
        // Absorbés (rouge)
        latex.SetTextColor(colAbsorbes);
        double yAbsorbes = (aAbsorbes[i] > 0) ? aAbsorbes[i] * 1.5 : 0.8;
        latex.DrawLatex(xAbsorbes - 0.1, yAbsorbes, Form("%.0f", aAbsorbes[i]));
    }
    latex.SetTextColor(kBlack);
    
    // Légende
    TLegend* leg5 = new TLegend(0.70, 0.75, 0.92, 0.88);
    leg5->AddEntry(hEmisParRaie, "Emis", "f");
    leg5->AddEntry(hEntresParRaie, "Entres eau", "f");
    leg5->AddEntry(hAbsorbesParRaie, "Absorbes", "f");
    leg5->Draw();
    
    c3->SaveAs("absorption_visualisation_recommandee.png");
    
    // ═════════════════════════════════════════════════════════════════════════
    // RÉSUMÉ FINAL
    // ═════════════════════════════════════════════════════════════════════════
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "  FICHIERS GENERES:" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << "  1. absorption_probleme_visualisation.png  - Comparaison des methodes" << std::endl;
    std::cout << "  2. absorption_explication_compton.png     - Explication effet Compton" << std::endl;
    std::cout << "  3. absorption_visualisation_recommandee.png - Graphiques recommandes" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    std::cout << "\n  CONCLUSION:" << std::endl;
    std::cout << "  -----------" << std::endl;
    std::cout << "  Les histogrammes en energie (hGammaEnteringWater) sont POLLUES" << std::endl;
    std::cout << "  par les gammas Compton diffuses dans l'air." << std::endl;
    std::cout << "  -> Utiliser les statistiques par RAIE (TTree gamma_lines)" << std::endl;
    std::cout << "     qui identifient chaque gamma par son energie INITIALE.\n" << std::endl;
    
    file->Close();
}
