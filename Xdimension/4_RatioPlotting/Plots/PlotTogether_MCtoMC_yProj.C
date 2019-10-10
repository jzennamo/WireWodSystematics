void PlotTogether_MCtoMC_yProj()
{
  gStyle->SetOptStat(0);
  TFile data_file("../../3_TreePlotting/out_hists_mc.root","read");
  TFile mc_file("../../3_TreePlotting/out_hists_mc_test.root","read");
  TH2D *dataHist;
  TH2D *mcHist;
  TCanvas *c_data = (TCanvas*)data_file.Get("hist_trk_x_hit_Q_plane_2_colz");
  TCanvas *c_mc   = (TCanvas*)mc_file.Get("hist_trk_x_hit_Q_plane_2_colz");
  dataHist = (TH2D*)c_data->GetPrimitive("hist_trk_x_hit_Q_plane_2");
  mcHist = (TH2D*)c_mc->GetPrimitive("hist_trk_x_hit_Q_plane_2");

  

  TH1D* mc_1D = (TH1D*)mcHist->ProjectionY();
  TH1D* data_1D = (TH1D*)dataHist->ProjectionY();
  
  TCanvas* c = new TCanvas("c");
  mc_1D->Draw("");
  data_1D->Draw("Same");

  TLegend* leg = new TLegend(0.4, 0.7, 0.85,0.85);
  leg->AddEntry(dataHist,"MC", "l");
  leg->AddEntry(mcHist,"MC - Wire Mod", "l");
  leg->Draw("same");
}
