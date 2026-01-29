// =============================================================================
// SCRIPT D'ANALYSE - PUITS COURONNE (Am-241)
// =============================================================================
// 
// Usage: root -l analyse_Dose.C
//    ou: root -l 'analyse_Dose.C("autre_fichier.root")'
//
// Adapté pour la simulation Am-241 (raies X et gamma 11-125 keV)
// =============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TString.h>
#include <iostream>
#include <iomanip>

// Couleurs pour les anneaux
const int ringColors[] = {kRed, kOrange+1, kGreen+2, kBlue, kViolet};

// Noms des raies gamma Am-241 (12 raies)
const char* gammaLineNames[] = {
    "X_{Ll} 11.9", "X_{L#alpha} 13.9", "X_{L#beta} 17.0", "X_{L#gamma} 20.8",
    "#gamma 26.3", "#gamma 33.2", "#gamma 43.4", "#gamma 55.6",
    "#gamma 59.5", "#gamma 99.0", "#gamma 103.0", "#gamma 125.3"
};
const int nGammaLines = 12;

// Énergies des raies Am-241 pour annotation
const double gammaEnergies[] = {11.89, 13.9, 17.0, 20.8, 26.3, 33.2, 43.4, 55.6, 59.5, 99.0, 103.0, 125.3};

void analyse_Dose(const char* filename = "output.root")
{
    // Configuration du style
    gStyle->SetOptStat(0);
    gStyle->SetHistLineWidth(3);
    gStyle->SetTitleFont(62, "");
    gStyle->SetTitleFontSize(0.06);
    
    // Nom de base pour les sorties
    TString baseName = filename;
    baseName.ReplaceAll(".root", "");
    
    // ==========================================================================
    // OUVERTURE DU FICHIER
    // ==========================================================================
    
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "ERREUR: Impossible d'ouvrir " << filename << std::endl;
        std::cerr << "Verifiez que le fichier existe avec: ls -la *.root" << std::endl;
        return;
    }
    
    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "     ANALYSE DU FICHIER: " << filename << "\n";
    std::cout << "     Source Am-241 (raies X/gamma 11-125 keV)\n";
    std::cout << "================================================================\n\n";
    
    // Lister le contenu du fichier
    std::cout << "Contenu du fichier:\n";
    file->ls();
    std::cout << "\n";
    
    // ==========================================================================
    // 1. SPECTRE DES GAMMAS ÉMIS
    // ==========================================================================
    
    TH1D* hGammaEmitted = (TH1D*)file->Get("hGammaEmitted");
    if (hGammaEmitted && hGammaEmitted->GetEntries() > 0) {
        TCanvas* c_spectrum = new TCanvas("c_spectrum", "Spectre gamma Am-241", 1200, 700);
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.05);
        
        // Ne pas rebinner pour Am-241 (raies fines)
        // hGammaEmitted->Rebin(2);

        hGammaEmitted->SetLineColor(kRed+2);
        hGammaEmitted->SetFillColor(kRed);
        hGammaEmitted->SetFillStyle(1001);
        hGammaEmitted->SetLineWidth(1);
        hGammaEmitted->GetXaxis()->SetRangeUser(0, 150);  // Am-241: 0-150 keV
        hGammaEmitted->GetXaxis()->SetTitle("Energie [keV]");
        hGammaEmitted->GetYaxis()->SetTitle("Counts");
        hGammaEmitted->SetTitle("Spectre gamma Am-241 emis");
        hGammaEmitted->Draw("HIST");

        // Annoter les raies principales Am-241
        TLatex* latex = new TLatex();
        latex->SetTextSize(0.03);
        latex->SetTextColor(kBlue+2);
        
        // Raies principales à annoter
        double energiesToLabel[] = {13.9, 17.0, 26.3, 59.5};
        const char* labelNames[] = {"X_{L#alpha}", "X_{L#beta}", "#gamma 26", "#gamma 59.5"};
        int nLabels = 4;
        
        for (int i = 0; i < nLabels; i++) {
            int bin = hGammaEmitted->FindBin(energiesToLabel[i]);
            double y = hGammaEmitted->GetBinContent(bin);
            if (y > 0) {
                latex->DrawLatex(energiesToLabel[i] + 2, y*1.8, labelNames[i]);
            }
        }
        
        c_spectrum->Update();
        c_spectrum->SaveAs(baseName + "_spectre_gamma_emis.png");
        std::cout << "=> Sauvegarde: " << baseName << "_spectre_gamma_emis.png\n";
        std::cout << "Spectre gamma emis: " << hGammaEmitted->GetEntries() << " entries\n\n";
    } else {
        std::cout << "ATTENTION: Histogramme hGammaEmitted vide ou non trouve!\n\n";
    }
    
    // ==========================================================================
    // 2. SPECTRE DES GAMMAS ENTRANT DANS L'EAU
    // ==========================================================================
    
    TH1D* hGammaWater = (TH1D*)file->Get("hGammaEnteringWater");
    if (hGammaWater && hGammaWater->GetEntries() > 0) {
        TCanvas* c_water = new TCanvas("c_water", "Spectre gamma entrant eau", 1200, 700);
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        
        hGammaWater->SetLineColor(kBlue+1);
        hGammaWater->SetLineWidth(2);
        hGammaWater->SetMarkerColor(kBlue+1);
        hGammaWater->SetMarkerStyle(20);
        hGammaWater->SetMarkerSize(0.8);
        hGammaWater->GetXaxis()->SetRangeUser(0, 150);  // Am-241: 0-150 keV
        hGammaWater->GetXaxis()->SetTitle("Energie [keV]");
        hGammaWater->GetYaxis()->SetTitle("Counts");
        hGammaWater->SetTitle("Spectre gamma Am-241 entrant dans l'eau");
        hGammaWater->Draw("P");
        
        c_water->Update();
        c_water->SaveAs(baseName + "_spectre_gamma_eau.png");
        std::cout << "=> Sauvegarde: " << baseName << "_spectre_gamma_eau.png\n";
        std::cout << "Spectre gamma eau: " << hGammaWater->GetEntries() << " entries\n\n";
    }
    
    // ==========================================================================
    // 3. HISTOGRAMMES DE DOSE PAR ANNEAU (nGy par événement)
    // ==========================================================================
    
    TCanvas* c_dose = new TCanvas("c_dose", "Dose par anneau", 1200, 800);
    c_dose->Divide(3, 2);
    
    TH1D* h_dose_ring[5];
    TH1D* h_dose_total = nullptr;
    
    bool hasHistos = false;
    
    for (int i = 0; i < 5; i++) {
        TString hname = Form("h_dose_ring%d", i);
        h_dose_ring[i] = (TH1D*)file->Get(hname);
        
        if (h_dose_ring[i] && h_dose_ring[i]->GetEntries() > 0) {
            hasHistos = true;
            c_dose->cd(i + 1);
            gPad->SetLogy();
            
            h_dose_ring[i]->SetLineColor(ringColors[i]);
            h_dose_ring[i]->SetMarkerColor(ringColors[i]);
            h_dose_ring[i]->SetMarkerStyle(21);
            h_dose_ring[i]->SetMarkerSize(0.6);
            h_dose_ring[i]->GetXaxis()->SetRangeUser(0, 0.15);  // Adapté Am-241
            h_dose_ring[i]->GetXaxis()->SetTitle("Dose [nGy]");
            h_dose_ring[i]->GetYaxis()->SetTitle("Counts");
            h_dose_ring[i]->GetXaxis()->SetTitleSize(0.05);
            h_dose_ring[i]->GetYaxis()->SetTitleSize(0.05);
            h_dose_ring[i]->GetXaxis()->SetLabelSize(0.05);
            h_dose_ring[i]->GetYaxis()->SetLabelSize(0.05);
            h_dose_ring[i]->SetTitle(Form("Anneau %d (r=%d-%d mm)", i, i*5, (i+1)*5));
            h_dose_ring[i]->Draw("P");
            
            gPad->Update();
            
            std::cout << "Anneau " << i << ": Entries=" << h_dose_ring[i]->GetEntries()
                      << ", Mean=" << std::scientific << h_dose_ring[i]->GetMean() << " nGy\n";
        }
    }

    // Dose totale
    h_dose_total = (TH1D*)file->Get("h_dose_total");
    if (h_dose_total && h_dose_total->GetEntries() > 0) {
        c_dose->cd(6);
        gPad->SetLogy();
        h_dose_total->SetLineColor(kBlack);
        h_dose_total->SetLineWidth(3);
        h_dose_total->SetMarkerColor(kBlack);
        h_dose_total->SetMarkerStyle(21);
        h_dose_total->SetMarkerSize(0.6);
        h_dose_total->GetXaxis()->SetRangeUser(0, 0.15);  // Adapté Am-241
        h_dose_total->GetXaxis()->SetTitle("Dose totale [nGy]");
        h_dose_total->GetYaxis()->SetTitle("Counts");
        h_dose_total->GetXaxis()->SetTitleSize(0.05);
        h_dose_total->GetYaxis()->SetTitleSize(0.05);
        h_dose_total->GetXaxis()->SetLabelSize(0.05);
        h_dose_total->GetYaxis()->SetLabelSize(0.05);
        h_dose_total->SetTitle("TOTAL (5 anneaux)");
        h_dose_total->Draw("P");

        gPad->Update();
        
        std::cout << "TOTAL: Entries=" << h_dose_total->GetEntries()
                  << ", Mean=" << std::scientific << h_dose_total->GetMean() << " nGy\n\n";
    }
    
    if (hasHistos) {
        c_dose->Update();
        c_dose->SaveAs(baseName + "_dose_par_anneau.png");
        std::cout << "=> Sauvegarde: " << baseName << "_dose_par_anneau.png\n\n";
    } else {
        std::cout << "ATTENTION: Aucun histogramme h_dose_ringX trouve avec des donnees!\n";
        std::cout << "Verifiez que la simulation a produit des depots d'energie.\n\n";
    }
    
    // ==========================================================================
    // 4. ENERGIE DEPOSEE PAR STEP (hEdepRing*)
    // ==========================================================================
    
    TCanvas* c_edep = new TCanvas("c_edep", "Energie deposee par step", 1200, 800);
    c_edep->Divide(3, 2);
    
    double statEdepX1 = 0.60, statEdepX2 = 0.98;
    double statEdepY1 = 0.82, statEdepY2 = 0.95;
    double statEdepTextSize = 0.05;
    
    bool hasEdep = false;
    for (int i = 0; i < 5; i++) {
        TString hname = Form("hEdepRing%d", i);
        TH1D* h = (TH1D*)file->Get(hname);
        
        if (h && h->GetEntries() > 0) {
            hasEdep = true;
            c_edep->cd(i + 1);
            gPad->SetLogy();
            
            h->SetLineColor(ringColors[i]);
            h->SetLineWidth(2);
            h->GetXaxis()->SetRangeUser(0, 80);  // Adapté Am-241 (max ~60 keV)
            h->GetXaxis()->SetTitle("Edep [keV]");
            h->GetYaxis()->SetTitle("Counts");
            h->SetTitle(Form("Edep/step Anneau %d", i));
            h->Draw();
            
            gPad->Update();
            TPaveStats* stats = (TPaveStats*)h->FindObject("stats");
            if (stats) {
                stats->SetX1NDC(statEdepX1);
                stats->SetX2NDC(statEdepX2);
                stats->SetY1NDC(statEdepY1);
                stats->SetY2NDC(statEdepY2);
                stats->SetTextSize(statEdepTextSize);
                stats->SetTextFont(62);
            }
            gPad->Modified();
            gPad->Update();
            
            std::cout << "hEdepRing" << i << ": " << h->GetEntries() << " entries, Mean=" 
                      << h->GetMean() << " keV\n";
        }
    }
    
    // Edep total dans l'eau
    TH1D* hEdepWater = (TH1D*)file->Get("hEdepWater");
    if (hEdepWater && hEdepWater->GetEntries() > 0) {
        c_edep->cd(6);
        gPad->SetLogy();
        hEdepWater->SetLineColor(kBlack);
        hEdepWater->SetLineWidth(2);
        hEdepWater->GetXaxis()->SetRangeUser(0, 80);  // Adapté Am-241
        hEdepWater->SetTitle("Edep total eau");
        hEdepWater->Draw();
        
        gPad->Update();
        TPaveStats* stats = (TPaveStats*)hEdepWater->FindObject("stats");
        if (stats) {
            stats->SetX1NDC(statEdepX1);
            stats->SetX2NDC(statEdepX2);
            stats->SetY1NDC(statEdepY1);
            stats->SetY2NDC(statEdepY2);
            stats->SetTextSize(statEdepTextSize);
            stats->SetTextFont(62);
        }
        gPad->Modified();
        gPad->Update();
        
        std::cout << "hEdepWater: " << hEdepWater->GetEntries() << " entries\n";
    }
    
    if (hasEdep) {
        c_edep->Update();
        c_edep->SaveAs(baseName + "_edep_par_step.png");
        std::cout << "=> Sauvegarde: " << baseName << "_edep_par_step.png\n\n";
    }
    
    // ==========================================================================
    // 5. TAUX D'ABSORPTION DANS L'EAU PAR RAIE GAMMA
    // ==========================================================================
    
    gStyle->SetOptStat(0);
    
    TTree* gamma_lines = (TTree*)file->Get("gamma_lines");
    if (gamma_lines) {
        std::cout << "\nAnalyse du TTree gamma_lines...\n";
        
        Int_t lineIndex;
        Double_t energy_keV;
        Int_t emitted, enteredWater, absorbedWater;
        Double_t waterAbsRate, waterEntryRate;
        
        gamma_lines->SetBranchAddress("lineIndex", &lineIndex);
        gamma_lines->SetBranchAddress("energy_keV", &energy_keV);
        gamma_lines->SetBranchAddress("emitted", &emitted);
        gamma_lines->SetBranchAddress("enteredWater", &enteredWater);
        gamma_lines->SetBranchAddress("absorbedWater", &absorbedWater);
        gamma_lines->SetBranchAddress("waterAbsRate", &waterAbsRate);
        gamma_lines->SetBranchAddress("waterEntryRate", &waterEntryRate);
        
        TCanvas* c_abs = new TCanvas("c_abs", "Taux d'absorption par raie", 1200, 700);
        c_abs->SetBottomMargin(0.20);
        c_abs->SetLeftMargin(0.12);
        c_abs->SetRightMargin(0.05);
        gPad->SetLogy();
        gPad->SetGridy();
        
        TH1D* h_abs_water = new TH1D("h_abs_water", 
            "Taux d'absorption dans l'eau (1 mm anneaux) par raie Am-241;Raie gamma;Taux d'absorption (%)", 
            nGammaLines, 0, nGammaLines);
        
        Long64_t nEntries = gamma_lines->GetEntries();
        std::cout << "Nombre d'entrees: " << nEntries << "\n";
        
        std::cout << "\n┌────────┬────────────┬───────────┬─────────────┬────────────┬──────────────┐\n";
        std::cout << "│ Index  │ Energie    │   Emis    │ Entré eau   │ Absorbé    │ Taux abs (%) │\n";
        std::cout << "├────────┼────────────┼───────────┼─────────────┼────────────┼──────────────┤\n";
        
        double maxRate = 0;
        for (Long64_t i = 0; i < nEntries && i < nGammaLines; i++) {
            gamma_lines->GetEntry(i);
            h_abs_water->SetBinContent(lineIndex + 1, waterAbsRate);
            h_abs_water->GetXaxis()->SetBinLabel(lineIndex + 1, gammaLineNames[lineIndex]);
            
            if (waterAbsRate > maxRate) maxRate = waterAbsRate;
            
            std::cout << "│   " << std::setw(2) << lineIndex << "   │ "
                      << std::setw(8) << std::fixed << std::setprecision(1) << energy_keV << " keV │"
                      << std::setw(10) << emitted << " │"
                      << std::setw(12) << enteredWater << " │"
                      << std::setw(11) << absorbedWater << " │"
                      << std::setw(12) << std::setprecision(2) << waterAbsRate << " │\n";
        }
        std::cout << "└────────┴────────────┴───────────┴─────────────┴────────────┴──────────────┘\n\n";
        
        h_abs_water->SetFillColor(kOrange+1);
        h_abs_water->SetLineColor(kOrange+3);
        h_abs_water->SetLineWidth(2);
        h_abs_water->SetBarWidth(0.85);
        h_abs_water->SetBarOffset(0.075);
        h_abs_water->SetMinimum(0.01);  // Minimum pour échelle log
        h_abs_water->SetMaximum(100);   // Max 100%
        
        h_abs_water->GetXaxis()->SetLabelSize(0.045);
        h_abs_water->GetXaxis()->SetLabelOffset(0.02);
        h_abs_water->GetXaxis()->SetTitleOffset(2.5);
        h_abs_water->GetYaxis()->SetLabelSize(0.045);
        h_abs_water->GetYaxis()->SetTitleSize(0.05);
        h_abs_water->GetYaxis()->SetTitleOffset(1.0);
        
        h_abs_water->Draw("bar");
        
        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.035);
        latex->DrawLatex(0.15, 0.87, "Source Am-241 (42 kBq) - Eau 1 mm (anneaux)");
        latex->DrawLatex(0.15, 0.82, Form("Absorption max: %.1f%% (raie X_{Ll} 11.9 keV)", maxRate));
        
        // Ligne horizontale à 50%
        TLine* line50 = new TLine(0, 50, nGammaLines, 50);
        line50->SetLineColor(kRed);
        line50->SetLineStyle(2);
        line50->SetLineWidth(2);
        line50->Draw();
        
        c_abs->Update();
        c_abs->SaveAs(baseName + "_taux_absorption_eau.png");
        std::cout << "=> Sauvegarde: " << baseName << "_taux_absorption_eau.png\n\n";
    } else {
        std::cout << "TTree gamma_lines non trouve dans le fichier\n";
        std::cout << "Le fichier ROOT a ete genere avec une ancienne version.\n\n";
    }
    
    // ==========================================================================
    // 6. CARTES 2D DE DEPOT D'ENERGIE
    // ==========================================================================
    
    gStyle->SetOptStat(0);
    
    TH2D* hEdepXY = (TH2D*)file->Get("hEdepXY");
    TH2D* hEdepRZ = (TH2D*)file->Get("hEdepRZ");
    
    if ((hEdepXY && hEdepXY->GetEntries() > 0) || (hEdepRZ && hEdepRZ->GetEntries() > 0)) {
        TCanvas* c_2d = new TCanvas("c_2d", "Cartes 2D", 1400, 600);
        c_2d->Divide(2, 1);
        
        if (hEdepXY && hEdepXY->GetEntries() > 0) {
            c_2d->cd(1);
            gPad->SetLogz();
            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.12);
            hEdepXY->SetTitle("Depot d'energie XY (Am-241);X [mm];Y [mm]");
            hEdepXY->Draw("COLZ");
            
            gPad->Update();
            std::cout << "hEdepXY: " << hEdepXY->GetEntries() << " entries\n";
        }
        
        if (hEdepRZ && hEdepRZ->GetEntries() > 0) {
            c_2d->cd(2);
            gPad->SetLogz();
            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.12);
            hEdepRZ->SetTitle("Depot d'energie RZ (Am-241);R [mm];Z [mm]");
            hEdepRZ->Draw("COLZ");
            
            std::cout << "hEdepRZ: " << hEdepRZ->GetEntries() << " entries\n";
        }
        
        c_2d->Update();
        c_2d->SaveAs(baseName + "_cartes_2d.png");
        std::cout << "=> Sauvegarde: " << baseName << "_cartes_2d.png\n\n";
    }
    
    // ==========================================================================
    // 7. SPECTRE DES ELECTRONS SECONDAIRES
    // ==========================================================================
    
    TH1D* hElectron = (TH1D*)file->Get("hElectronSpectrum");
    if (hElectron && hElectron->GetEntries() > 0) {
        TCanvas* c_elec = new TCanvas("c_elec", "Electrons secondaires", 1000, 700);
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        
        hElectron->SetLineColor(kGreen+2);
        hElectron->SetLineWidth(2);
        hElectron->SetMarkerColor(kGreen+2);
        hElectron->SetMarkerStyle(20);
        hElectron->SetMarkerSize(0.6);
        hElectron->GetXaxis()->SetRangeUser(0, 80);  // Adapté Am-241
        hElectron->SetTitle("Spectre des electrons secondaires dans l'eau (Am-241)");
        hElectron->Draw("P");
        
        c_elec->Update();
        c_elec->SaveAs(baseName + "_spectre_electrons.png");
        std::cout << "=> Sauvegarde: " << baseName << "_spectre_electrons.png\n";
        std::cout << "Electrons: " << hElectron->GetEntries() << " entries, Mean=" 
                  << hElectron->GetMean() << " keV\n\n";
    }
    
    // ==========================================================================
    // 8. COMPARAISON SPECTRES ÉMIS vs ENTRÉ vs ABSORBÉ
    // ==========================================================================
    
    if (hGammaEmitted && hGammaWater) {
        TCanvas* c_comp = new TCanvas("c_comp", "Comparaison spectres", 1200, 700);
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.05);
        
        // Cloner pour ne pas modifier les originaux
        TH1D* hEmis = (TH1D*)hGammaEmitted->Clone("hEmis_clone");
        TH1D* hEntres = (TH1D*)hGammaWater->Clone("hEntres_clone");
        
        hEmis->SetLineColor(kBlue+1);
        hEmis->SetLineWidth(2);
        hEmis->SetFillStyle(0);
        hEmis->GetXaxis()->SetRangeUser(0, 150);
        hEmis->SetTitle("Comparaison spectres Am-241: Emis vs Entres eau");
        hEmis->Draw("HIST");
        
        hEntres->SetLineColor(kGreen+2);
        hEntres->SetLineWidth(2);
        hEntres->SetFillStyle(0);
        hEntres->Draw("HIST SAME");
        
        TLegend* leg = new TLegend(0.55, 0.75, 0.90, 0.88);
        leg->SetBorderSize(1);
        leg->SetFillColor(kWhite);
        leg->AddEntry(hEmis, Form("Emis (N=%d)", (int)hEmis->GetEntries()), "l");
        leg->AddEntry(hEntres, Form("Entres eau (N=%d)", (int)hEntres->GetEntries()), "l");
        leg->Draw();
        
        c_comp->Update();
        c_comp->SaveAs(baseName + "_comparaison_spectres.png");
        std::cout << "=> Sauvegarde: " << baseName << "_comparaison_spectres.png\n\n";
    }
    
    // ==========================================================================
    // RESUME
    // ==========================================================================
    
    std::cout << "================================================================\n";
    std::cout << "                    ANALYSE TERMINEE                            \n";
    std::cout << "================================================================\n";
    std::cout << "  Fichiers generes:                                             \n";
    if (hGammaEmitted && hGammaEmitted->GetEntries() > 0)
        std::cout << "    - " << baseName << "_spectre_gamma_emis.png\n";
    if (hGammaWater && hGammaWater->GetEntries() > 0)
        std::cout << "    - " << baseName << "_spectre_gamma_eau.png\n";
    if (hasHistos)
        std::cout << "    - " << baseName << "_dose_par_anneau.png\n";
    if (hasEdep)
        std::cout << "    - " << baseName << "_edep_par_step.png\n";
    if (gamma_lines)
        std::cout << "    - " << baseName << "_taux_absorption_eau.png\n";
    if ((hEdepXY && hEdepXY->GetEntries() > 0) || (hEdepRZ && hEdepRZ->GetEntries() > 0))
        std::cout << "    - " << baseName << "_cartes_2d.png\n";
    if (hElectron && hElectron->GetEntries() > 0)
        std::cout << "    - " << baseName << "_spectre_electrons.png\n";
    if (hGammaEmitted && hGammaWater)
        std::cout << "    - " << baseName << "_comparaison_spectres.png\n";
    std::cout << "================================================================\n\n";
    
    // Ne pas fermer le fichier pour permettre l'exploration interactive
    // file->Close();
}
