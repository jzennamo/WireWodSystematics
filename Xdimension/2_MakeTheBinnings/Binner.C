void Binner(){

  TFile* f = new TFile("../1_TreeMaking/outfile_mc_DLmod.root");
  TTree* t = (TTree*)f->Get("hit_v_track");

  std::vector<string> trk_kin  = {"trk_x"};//, "trk_ThetaXZ"}; //,"trk_ThetaYZ","trk_phi"*/}; 
  std::vector<double> trk_low  = {1.6 };//, -3.2}; //, -3.2, -3.2*/};
  std::vector<double> trk_high = {250.};//, 3.2}; //,  3.2,  3.2*/};
  std::vector<std::vector<std::vector< int > > > bins; 
  std::vector<std::vector<std::vector <double> > > edges;
  std::vector<std::vector<TH2D* > > hist;
  bins.resize(3);
  edges.resize(3);
  hist.resize(3);
  for(int pl = 0; pl < 3; pl++){
    bins[pl].resize(trk_kin.size());
    edges[pl].resize(trk_kin.size());
    hist[pl].resize(trk_kin.size());
    for(int trk = 0; trk < trk_kin.size(); trk++){
      
      hist[pl][trk] = new TH2D(Form("hist_%s_%d", trk_kin[trk].c_str(), pl),"", 1000, trk_low[trk], trk_high[trk], 1000, -1 , 1);
      
      t->Project(Form("hist_%s_%d", trk_kin[trk].c_str(), pl),Form("trk_cosTheta:%s",trk_kin[trk].c_str()),Form("(1/Nhits)*(hit_plane == %d && abs(hit_sigma - 5.333333) > 0.0001)",pl));
 
      int leadbin = 1;
      double sum = 0;

      bins[pl][trk].push_back(1);
               
      for(int i = 1; i <= hist[pl][trk]->GetNbinsX(); i++){
	
	if(hist[pl][trk]->Integral(leadbin, i, 0, 1000) >= 1000){
	  
	  leadbin = i;
	  bins[pl][trk].push_back(i);
	}    

	
      }
    }
  }

  for(int pl = 0; pl < 3; pl++){
    for(int trk = 0; trk < trk_kin.size(); trk++){
           
      std::cout << "Variable : " << trk_kin[trk].c_str() << " , Plane : " << pl << std::endl; 
      for(int i = 0; i < bins[pl][trk].size(); i++){
     
	std::cout << hist[pl][trk]->GetXaxis()->GetBinLowEdge(bins[pl][trk][i]) << ", "; 
	edges[pl][trk].push_back(hist[pl][trk]->GetXaxis()->GetBinLowEdge(bins[pl][trk][i]));
      }
      std::cout << "\n" << std::endl;
    }
  }
  
    
  std::vector< std::vector < std::vector<  TH1D* > > > all_da_Plots;
  std::vector< std::vector < TLegend* > > leg;
  all_da_Plots.resize(trk_kin.size());
  leg.resize(trk_kin.size());
  for( int i = 0; i < trk_kin.size(); i++){
    all_da_Plots[i].resize(3);
    leg[i].resize(3);
    for(int pl = 0; pl < 3; pl++){
      leg[i][pl] = new TLegend(0.4,0.5,0.8,0.8);
      leg[i][pl]->SetNColumns(2);
      all_da_Plots[i][pl].resize(edges[pl][i].size()-1);
    }
  }
  
  std::cout << "Filling the plots... " << std::endl;
  
  for( int k = 0; k < trk_kin.size(); k++){
    for(int pl = 0; pl < 3; pl++){
      for(int i = 0; i < edges[pl][k].size()-1; i++){

	std::cout << " k " << k << "  pl " << pl << " edge  " << i << std::endl;
	
	all_da_Plots[k][pl][i] = new TH1D(Form("hist_%s_%d_bin_%d", trk_kin[k].c_str(), pl,i),"",30,-1,1);
	t->Project(Form("hist_%s_%d_bin_%d", trk_kin[k].c_str(), pl,i),
		   "trk_cosTheta",
		   Form("(1/Nhits)*(hit_plane == %d && abs(hit_sigma - 5.333333) > 0.0001 && trk_x > %f && trk_x < %f)", pl, edges[pl][k][i],edges[pl][k][i+1]));
	
	leg[k][pl]->AddEntry(all_da_Plots[k][pl][i],Form("[%4.2f cm, %4.2f cm]",edges[pl][k][i],edges[pl][k][i+1]),"l");
      }
    }
  }

  std::cout << "Plotting the plots..." << std::endl;

  std::vector<std::vector< TCanvas* > > c;
  c.resize(trk_kin.size());
  for( int k = 0; k < trk_kin.size(); k++){
    c[k].resize(3);
    for(int pl = 0; pl < 3; pl++){      
      std::cout << "make the canvas..." << std::endl;
      c[k][pl] = new TCanvas(Form("c_%s_%d", trk_kin[k].c_str(), pl));
      c[k][pl]->cd();
      all_da_Plots[k][pl][0]->SetLineColor(1);
      all_da_Plots[k][pl][0]->SetLineWidth(2);
      std::cout << "...draw the plots..." << std::endl;
      all_da_Plots[k][pl][0]->DrawNormalized("hist",1);
      for(int i = 1; i < edges[pl][k].size()-1; i++){	
	if(i <= 8){
	  all_da_Plots[k][pl][i]->SetLineColor(i+1);
	}
	else{
	  all_da_Plots[k][pl][i]->SetLineColor(40+(i+1)-8);
	}
	all_da_Plots[k][pl][i]->SetLineWidth(2);
	all_da_Plots[k][pl][i]->DrawNormalized("histsame",1);
      }
      std::cout << "...and legend." << std::endl;
      leg[k][pl]->Draw("same");
    }
  }
}
