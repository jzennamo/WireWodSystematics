void PlotsOfData(){

  TFile* f = new TFile("../1_TreeMaking/outfile_mc_DLmod.root");
  TTree* t = (TTree*)f->Get("hit_v_track");

  TFile* fd = new TFile("../1_TreeMaking/outfile_mc_DLmod.root");
  TTree* td = (TTree*)fd->Get("hit_v_track");

  TH2D* hist = new TH2D("hist","", 1000, 0, 256, 1000, -1 , 1);

  t->Project("hist","trk_cosTheta:trk_x","(1/Nhits)*(hit_plane == 0 && abs(hit_sigma - 5.333333) > 0.0001 && max(trk_StartX,trk_EndX) < 256 && min(trk_StartX,trk_EndX) > 0.3)");
 
  int leadbin = 1;
  double sum = 0;
  std::vector< int > bins; 
  bins.push_back(1);
  for(int i = 1; i <= hist->GetNbinsX(); i++){
    
    if(hist->Integral(leadbin, i, 0, 1000) >= 2500){
      std::cout << " from " << leadbin << " to " << i << " has " << hist->Integral(leadbin, i, 0, 1000) << std::endl; 
      leadbin = i;
      bins.push_back(i);
    }    
  }

  std::cout << " Bins " << std::endl; 
  std::vector <double> edges;
  for(int i = 0; i < bins.size(); i++){
    
    std::cout << hist->GetXaxis()->GetBinLowEdge(bins[i]) << std::endl;
    edges.push_back(hist->GetXaxis()->GetBinLowEdge(bins[i]));
  }
  

  std::vector< TH1D* > all_da_Plots;
  all_da_Plots.resize(edges.size()-1);

  TLegend* leg = new TLegend(0.4,0.5,0.8,0.8);
  leg->SetNColumns(2);
  for(int i = 0; i < edges.size()-1; i++){
    all_da_Plots[i] = new TH1D(Form("bin_%d",i),"",30,-1,1);
    td->Project(Form("bin_%d",i),
		"trk_cosTheta",
		Form("(1/Nhits)*(hit_plane == 0 && max(trk_StartX,trk_EndX) < 256 && min(trk_StartX,trk_EndX) > 0.3 && abs(hit_sigma - 5.333333) > 0.0001 && trk_x > %f && trk_x < %f)",edges[i],edges[i+1]));
    
    leg->AddEntry(all_da_Plots[i],Form("[%4.2f cm, %4.2f cm]",edges[i],edges[i+1]),"l");
  }
  
  TCanvas* c = new TCanvas("c");
  all_da_Plots[0]->SetLineColor(1);
  all_da_Plots[0]->SetLineWidth(2);
  all_da_Plots[0]->DrawNormalized("hist",1);
  for(int i = 1; i < edges.size()-1; i++){
    if(i <= 8){
      all_da_Plots[i]->SetLineColor(i+1);
    }
    else{
      all_da_Plots[i]->SetLineColor(40+(i+1)-8);
    }
    all_da_Plots[i]->SetLineWidth(2);
    all_da_Plots[i]->DrawNormalized("histsame",1);
  }
  leg->Draw("same");
}


