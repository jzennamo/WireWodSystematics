void PlotTogether_MCtoMC()
{
  gStyle->SetOptStat(0);
  TFile data_file("../../3_TreePlotting/out_hists_mc.root","read");
  TFile mc_file("../../3_TreePlotting/out_hists_mc_wiremod.root","read");
  TH1D *dataHist;
  TH1D *mcHist;
  TCanvas *c_data = (TCanvas*)data_file.Get("hist_trk_x_hit_Q_plane_2_prof");
  TCanvas *c_mc   = (TCanvas*)mc_file.Get("hist_trk_x_hit_Q_plane_2_prof");
  dataHist = (TH1D*)c_data->GetPrimitive("histtm_trk_x_hit_Q_plane_2");
  mcHist = (TH1D*)c_mc->GetPrimitive("histtm_trk_x_hit_Q_plane_2");

  TCanvas* c = new TCanvas("c");
  //  mcHist->Divide(dataHist);
  mcHist->SetLineColor(kRed);
  mcHist->SetLineWidth(2);
  //  mcHist->GetYaxis()->SetTitle("Hit Charge [ADC*ticks]");
  // mcHist->GetXaxis()->SetTitle("Track x Position [cm]");
  mcHist->Draw();
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(2);
  dataHist->Draw("same");

  TLegend* leg = new TLegend(0.4, 0.7, 0.85,0.85);
  leg->AddEntry(dataHist,"MC", "l");
  leg->AddEntry(mcHist,"MC - New Trun Mean", "l");
  leg->Draw("same");
}
