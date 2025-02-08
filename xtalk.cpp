#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>


#define savePlot

int nRows = 336;
int nColumns = 432;

using namespace std;

void xtalk(){
   
   // ---- ---- ---- ---- Retrieving files and histograms ---- ---- ---- ----
   
   TFile *f_injtype1 = new TFile("runs/injType1.root");
   TFile *f_injtype5 = new TFile("runs/injType5.root");
   TFile *f_injtype6 = new TFile("runs/injType6.root");
   
   
   float alive_eff = 0.7;
   float coupled_eff = 0.3;
   float uncoupled_eff = 0.1;
   
   string baseDir = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_";
   string shortBaseDir = "D_B(0)_O(0)_H(0)_";
   string moduleID = "P1002";
   string cycle = "00";
   
   const std::filesystem::path path = "plots/xtalk/"+moduleID+"_cycles/broken_bumps_"+moduleID+"_cycle"+cycle+".txt";
   std::ofstream broken_list(path);
   if (!broken_list.is_open()) return;
   
   
   TCanvas *c_full_quad_confirmed = new TCanvas("c_full_quad_confirmed",("confirmed broken bumps for module "+moduleID).c_str(),1600,1600);
   c_full_quad_confirmed->Divide(2,2);
   int chipMap[4] = {1,4,3,2};
   
   
   TCanvas *c_full_confirmed_map = new TCanvas("c_full_confirmed_map",("confirmed broken bumps for module "+moduleID).c_str(),nColumns*6,nRows*6);
   c_full_confirmed_map->SetRightMargin(0.14);
   TH2F *h_full_confirmed_map = new TH2F("h_full_confirmed_map",("confirmed broken bumps for module "+moduleID+";;(connector side)                    ").c_str(),nColumns*2,0,nColumns*2,nRows*2,0,nRows*2);
   h_full_confirmed_map->GetXaxis()->SetLabelOffset(0.17);
   h_full_confirmed_map->GetYaxis()->SetLabelOffset(0.17);
   h_full_confirmed_map->GetYaxis()->SetTitleSize(0.05);
   h_full_confirmed_map->GetYaxis()->SetTitleOffset(0.7);
   
   for (int i=0; i<2*nRows; i++){
       for (int j=0; j<2*nColumns; j++){
         h_full_confirmed_map->SetBinContent(j+1,i+1,0.0001);
       }
   }
   
   
   int chipID[4] = {15,14,12,13};
   
   //for (int ch = 15; ch >= 12; ch--){
   for (int pos = 0; pos < 4; pos++){
     int ch = chipID[pos];
     //if (ch==14) continue;
     string chip = to_string(ch);
     string plotDir = "plots/xtalk/"+moduleID+"_cycles/"+cycle+"/chip"+chip+"/";
     
     TCanvas *c0_pixelalive1 = (TCanvas*) f_injtype1->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive1 = (TH2F*)c0_pixelalive1->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TCanvas *c0_pixelalive5 = (TCanvas*) f_injtype5->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive5 = (TH2F*)c0_pixelalive5->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TCanvas *c0_pixelalive6 = (TCanvas*) f_injtype6->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive6 = (TH2F*)c0_pixelalive6->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     
     int nColumns = h_pixelalive5->GetXaxis()->GetNbins();
     int nRows = h_pixelalive5->GetYaxis()->GetNbins();
     
     //cout << "nRows: " << nRows << "  nColumns: " << nColumns << endl;
     
     
     
     
     
     
     // ---- ---- ---- ---- analysis ---- ---- ---- ----
     
     vector<int> dead_row,dead_col,suspicious_row,suspicious_col,confirmed_row,confirmed_col;
     TH2F *h_suspicious2D = (TH2F*)h_pixelalive5->Clone("h_suspicious2D");
     h_suspicious2D->SetTitle((moduleID+": suspicious channels of chip "+chip).c_str());
     TH2F *h_confirmed2D = (TH2F*)h_pixelalive5->Clone("h_confirmed2D");
     h_confirmed2D->SetTitle((moduleID+": confirmed disconnected channels of chip "+chip).c_str());
     
     //TH2F *h_zoomed_confirmed2D = new TH2F("h_zoomed_confirmed2D", ("confirmed disconnected channels of chip "+chip+";Columns;Rows").c_str(), 4,0,4, 250,0,250);
     
     cout << "chip " << ch << ":" << endl;
     broken_list << "chip " << ch << ":" << endl;
     
     
     bool mirror = false;
     if (pos>1) mirror = true;
     int offsetX = nColumns*(pos%2);
     int offsetY = nRows*(1-int(pos/2));
     
     
     for (int i=0; i<nRows; i++){
       for (int j=0; j<nColumns; j++){
         int newRow = i;
         int newCol = j;
         if(mirror) newRow = nRows-i-1;
         else newCol = nColumns-j-1;
         bool detectable = true;
         if (h_pixelalive1->GetBinContent(j+1,i+1)<alive_eff){
           dead_row.push_back(i);
           dead_col.push_back(j);
           detectable = false;
         }
         if (detectable){
           if ((i==0 && j%2 == 0) || (i==nRows-1 && j%2 == 1)) detectable = false;
         }
         if (detectable && h_pixelalive1->GetBinContent(j+1,i+1)>alive_eff && h_pixelalive5->GetBinContent(j+1,i+1)<coupled_eff){
           suspicious_row.push_back(i);
           suspicious_col.push_back(j);
           h_suspicious2D->SetBinContent(j+1,i+1,1);
         }
         else h_suspicious2D->SetBinContent(j+1,i+1,0.0001);
         if (detectable && h_pixelalive1->GetBinContent(j+1,i+1)>alive_eff && h_pixelalive5->GetBinContent(j+1,i+1)<coupled_eff && h_pixelalive6->GetBinContent(j+1,i+1)<uncoupled_eff){
           confirmed_row.push_back(i);
           confirmed_col.push_back(j);
           //cout << "row: " << i << " column: " << j << endl;
           broken_list << "row: " << i << " column: " << j << endl;
           h_confirmed2D->SetBinContent(j+1,i+1,1);
           h_full_confirmed_map->SetBinContent(offsetX+newCol+1,offsetY+newRow+1,1);
           //if (i<250 && j<4) h_zoomed_confirmed2D->SetBinContent(j+1,i+1,1);
         }
         else{
           h_confirmed2D->SetBinContent(j+1,i+1,0.0001);
           //h_full_confirmed_map->SetBinContent(offsetX+newCol+1,offsetY+newRow+1,0.0001);
           //if (i<250 && j<4) h_zoomed_confirmed2D->SetBinContent(j+1,i+1,0.0001);
         }
       }
     }
     
     cout << "    dead:        " << dead_row.size() << " (" << int(1000*dead_row.size()/(nRows*nColumns))/10. << "%)" << endl;
     cout << "    suspicious:  " << suspicious_row.size() << " (" << int(1000*suspicious_row.size()/(nRows*nColumns))/10. << "%)" << endl;
     cout << "    confirmed:   " << confirmed_row.size() << " (" << int(1000*confirmed_row.size()/(nRows*nColumns))/10. << "%)" << endl;
     
     broken_list << "    dead:        " << dead_row.size() << " (" << int(1000*dead_row.size()/(nRows*nColumns))/10. << "%)" << endl;
     broken_list << "    suspicious:  " << suspicious_row.size() << " (" << int(1000*suspicious_row.size()/(nRows*nColumns))/10. << "%)" << endl;
     broken_list << "    confirmed:   " << confirmed_row.size() << " (" << int(1000*confirmed_row.size()/(nRows*nColumns))/10. << "%)" << endl;
     
     
     
     
     
     
     
     
     
     
     // ---- ---- ---- ---- plotting ---- ---- ---- ----
     
     gStyle->SetOptStat(0);
     
     TCanvas *c_pixelalive1 = new TCanvas(("c_pixelalive1_"+chip).c_str(),("efficiency when injecting in same pixel of chip "+chip).c_str(),800,800);
     h_pixelalive1->Draw("colz");
     
     TCanvas *c_pixelalive5 = new TCanvas(("c_pixelalive5_"+chip).c_str(),("efficiency when injecting in coupled pixel of chip "+chip).c_str(),800,800);
     h_pixelalive5->Draw("colz");
     
     TCanvas *c_pixelalive6 = new TCanvas(("c_pixelalive6_"+chip).c_str(),("efficiency when injecting in uncoupled pixel of chip "+chip).c_str(),800,800);
     h_pixelalive6->Draw("colz");
     
     TCanvas *c_suspicious2D = new TCanvas(("c_suspicious2D_"+chip).c_str(),("suspicious channels of chip "+chip).c_str(),800,800);
     h_suspicious2D->Draw("colz");
     
     TCanvas *c_confirmed2D = new TCanvas(("c_confirmed2D_"+chip).c_str(),("confirmed disconnected channels of chip "+chip).c_str(),800,800);
     h_confirmed2D->Draw("colz");
     
     //TCanvas *c_zoomed_confirmed2D = new TCanvas(("c_zoomed_confirmed2D_"+chip).c_str(),("zoomed_confirmed disconnected channels of chip "+chip).c_str(),800,800);
     //h_zoomed_confirmed2D->Draw("colz");


#ifdef savePlot
     c_pixelalive1->SaveAs((plotDir+"pixelalive_"+chip+".png").c_str());
     c_pixelalive5->SaveAs((plotDir+"eff_coupled_"+chip+".png").c_str());
     c_pixelalive6->SaveAs((plotDir+"eff_uncoupled_"+chip+".png").c_str());
     c_suspicious2D->SaveAs((plotDir+"suspicious2D_"+chip+".png").c_str());
     c_confirmed2D->SaveAs((plotDir+"confirmed2D_"+chip+".png").c_str());
     //c_zoomed_confirmed2D->SaveAs((plotDir+"zoomed_confirmed2D_"+chip+".png").c_str());
#endif
     
     
     c_full_quad_confirmed->cd(chipMap[ch-12]);
     h_confirmed2D->Draw("colz");
     
     
   } //loop on chips
   
   //c_full_quad_confirmed->SaveAs(("plots/xtalk/"+moduleID+"_cycles/"+moduleID+"_cycle"+cycle+"_full_quad_confirmed.png").c_str());
   
   
   broken_list.close();
   
   
   c_full_confirmed_map->cd();
   //c_full_confirmed_map->SetLogz(1);
   h_full_confirmed_map->Draw("colz");
   c_full_confirmed_map->SaveAs(("plots/xtalk/"+moduleID+"_cycles/"+moduleID+"_cycle"+cycle+"_full_confirmed_map.png").c_str());
   
  
}

