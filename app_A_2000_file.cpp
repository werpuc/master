void app_A_2000_file(){

	/*
	events with max. 2 tracks per site
	- limitation to transverse momentum p_t
	- filling necessary and possibly unnecessary histograms
	- dividing them
	- all to file

	names are in polish for my thesis
	*/


	TChain *d3pd = new TChain("HeavyIonD3PD", "");
	HeavyIonD3PD hichain(d3pd, "data-slim-all.list");
	

	d3pd->SetBranchStatus("*", 0);
	d3pd->SetBranchStatus("trk_n", 1);
	d3pd->SetBranchStatus("trk_phi", 1);
	d3pd->SetBranchStatus("trk_eta", 1);
	d3pd->SetBranchStatus("trk_pt", 1);
	d3pd->SetBranchStatus("mb_E");

	gStyle->SetOptStat(0000000);

	int j; // number of events
	int all = 0; // all reconstructed tracks

	// temporary values	
	float eta;
	float phi;	

	float diff; 

	// to get data
	int tmp;
	float temp;
	
	int trkA[2]; // tracks numbers
	int trkC[2];
	int nA, nC;

	int nevents = 0; // number of tracks with min 2 tracks
	int neventsclear = 0; // without duplicates
	

	int mod; // help
	int m;
	
	

	float pt_cut = 2000; // cut 2 GeV

	TH1F *hist_A = new TH1F("hist_A", " ; #phi; liczba sladow", 90, -180, 180);
	TH1F *hist_C = new TH1F("hist_C", " ; #phi; liczba sladow", 90, -180, 180);

	TH1F *h_S[32];
	TH1F *h_B[32];

	for(int i=0; i<32; i++){
		char name[20];
		char names[50];
		sprintf(name, "h_S_%i", i);
		sprintf(names, " ; #phi; liczba sladow", i);
		h_S[i] = new TH1F(name, names, 90, -180, 180);
		sprintf(name, "h_B_%i", i);
		sprintf(names, " ; #phi; liczba sladow", i);
		h_B[i] = new TH1F(name, names, 90, -180, 180);
	}

	cout << "start" << endl;

	for(j=0; ; j++){
		int st = hichain.GetEntry(j);
		if(st<1) break;	
		// ************
		//getting data
		tmp = hichain.trk_n; 
		if(tmp!=0){
			temp = hichain.trk_eta->at(0);
			temp = hichain.trk_phi->at(0);
			temp = hichain.trk_pt->at(0);
			temp = hichain.mb_E->at(0);
		}		
		for(int i=0; i<2; i++){
			trkA[i] = -1;
			trkC[i] = -1;			
		}
		nA = nC = 0;
		for(int i=0; i<hichain.trk_n; i++){
			eta = hichain.trk_eta->at(i);
			if(eta>2.196){
				if(nA==0){
					trkA[0] = i;
					nA++;
				}
				else {
					trkA[1] = i;
					nA++;
				}
			}
			if(eta<-2.196){
				if(nC==0){
					trkC[0] = i;
					nC++;
				}
				else {
					trkC[1] = i;
					nC++;
				}
			}
		}
		if((trkA[0]+trkC[0])<0) continue;
		nevents++;
		// events with one track!
		if(trkA[1]!=-1){
			//diff = 57.296*fabs(hichain.trk_phi->at(trkA[0])-hichain.trk_phi->at(trkA[1]));
			//if(diff>180) diff-=180;
			//if(diff<60){
				trkA[0] = -1;
				trkA[1] = -1;
			//}
		}
		if(trkC[1]!=-1){
			//diff = 57.296*fabs(hichain.trk_phi->at(trkC[0])-hichain.trk_phi->at(trkC[1]));
			//if(diff>180) diff-=180;
			//if(diff<60){
				trkC[0] = -1;
				trkC[1] = -1;
			//}
		}
		if((trkA[0]+trkC[0])<0) continue;
		//neventsclear++;
		// filling histograms - depending if proper signal		
		for(int k=0; k<2; k++){
			if(trkA[k]!=-1){
				phi = 57.296*hichain.trk_phi->at(trkA[k]);
				if(hichain.trk_pt->at(trkA[k])>pt_cut){
					hist_A->Fill(phi);
					for(int l=16; l<32; l++){
						if(l%2!=0) continue;
						if(hichain.mb_E->at(l)>0.1) h_S[l]->Fill(phi);
						else h_B[l]->Fill(phi);
					}
				}
			}
			if(trkC[k]!=-1){
				phi = 57.296*hichain.trk_phi->at(trkC[k]);
				if(hichain.trk_pt->at(trkC[k])>pt_cut){
					hist_C->Fill(phi);
					for(int l=0; l<16; l++){
						if(l%2!=0) continue;
						if(hichain.mb_E->at(l)>0.1) h_S[l]->Fill(phi);
						else h_B[l]->Fill(phi);
					}
				}
			}
		}
	}

	cout << "the end" << endl;

	TGraphAsymmErrors *gr_SA[16];
	TGraphAsymmErrors *gr_BA[16];
	int temp_A_bins = hist_A->GetNbinsX();
	int temp_C_bins = hist_C->GetNbinsX();

	// calculating errors
	for(int i=0; i<8; i++){
		gr_SA[i] = new TGraphAsymmErrors(temp_C_bins);		
		gr_BA[i] = new TGraphAsymmErrors(temp_C_bins);
	}
	for(int i=8; i<16; i++){
		gr_SA[i] = new TGraphAsymmErrors(temp_A_bins);		
		gr_BA[i] = new TGraphAsymmErrors(temp_A_bins);
	}

	// dividing to get efficiency
	int temp_mod;
	for(int i=0; i<8; i++){
		temp_mod = 1+2*i;
		gr_SA[i]->BayesDivide(h_S[temp_mod], hist_C);
		gr_BA[i]->BayesDivide(h_B[temp_mod], hist_C);
	}
	for(int i=8; i<16; i++){
		temp_mod = 1+2*i;
		gr_SA[i]->BayesDivide(h_S[temp_mod], hist_A);
		gr_BA[i]->BayesDivide(h_B[temp_mod], hist_A);
	}

	// just for nice picture
	TLatex Tl;
    Tl.SetTextAlign(11);
    Tl.SetTextSize(0.03);

	TCanvas *c_mod_1_S_B = new TCanvas("c_mod_1_2000_ready_S_B", "c_mod_1_2000_ready_S_B", 1000, 600);
	c_mod_1_S_B->Divide(2, 1, 0.001, 0.001);
	c_mod_1_S_B->cd(1);
	h_S[1]->SetAxisRange(0., 400,"Y");
	h_S[1]->Draw();
	Tl.DrawLatex(-150, 380, "ATLAS work in progress");
    Tl.DrawLatex(-150, 360, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(-150, 340, "p_{t}>2000 GeV, |#eta|>2.196");
	c_mod_1_S_B->cd(2);
	h_B[1]->SetAxisRange(0., 400,"Y");
	h_B[1]->Draw();
	Tl.DrawLatex(-150, 380, "ATLAS work in progress");
    Tl.DrawLatex(-150, 360, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(-150, 340, "p_{t}>2000 GeV, |#eta|>2.196");
	c_mod_1_S_B->SaveAs("mod_1_SB_ready.C");
	c_mod_1_S_B->SaveAs("mod_1_SB_ready.gif");

	TCanvas *c_mod_17_S_B = new TCanvas("c_mod_17_2000_ready_S_B", "c_mod_17_2000_ready_S_B", 1000, 600);
	c_mod_17_S_B->Divide(2, 1, 0.001, 0.001);
	c_mod_17_S_B->cd(1);
	h_S[17]->SetAxisRange(0., 300,"Y");
	h_S[17]->Draw();
	Tl.DrawLatex(-150, 280, "ATLAS work in progress");
    Tl.DrawLatex(-150, 270, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(-150, 260, "p_{t}>2000 GeV, |#eta|>2.196");
	c_mod_17_S_B->cd(2);
	h_B[17]->SetAxisRange(0., 300,"Y");
	h_B[17]->Draw();
	Tl.DrawLatex(-150, 280, "ATLAS work in progress");
    Tl.DrawLatex(-150, 270, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(-150, 260, "p_{t}>2000 GeV, |#eta|>2.196");
	c_mod_17_S_B->SaveAs("mod_17_SB_ready.C");
	c_mod_17_S_B->SaveAs("mod_17_SB_ready.gif");

	TCanvas *c_norm_A = new TCanvas("c_norm_2000_A", "c_norm_2000_A", 1000, 600);
	hist_A->SetAxisRange(0., 300,"Y");
	hist_A->Draw("E");
	Tl.DrawLatex(-150, 280, "ATLAS work in progress");
    Tl.DrawLatex(-150, 270, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(-150, 260, "p_{t}>2000 GeV, |#eta|>2.196");
	c_norm_A->SaveAs("norm_2000_A.C");
	c_norm_A->SaveAs("norm_2000_A.gif");

	cout << "entries in A: " << hist_A->GetEntries() << endl;
	cout << "entries in C: " << hist_C->GetEntries() << endl;

	TCanvas *c_norm_C = new TCanvas("c_norm_2000_C", "c_norm_2000_C", 1000, 600);
	hist_C->SetAxisRange(0., 400,"Y");
	hist_C->Draw("E");
	Tl.DrawLatex(-150, 380, "ATLAS work in progress");
    Tl.DrawLatex(-150, 360, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(-150, 340, "p_{t}>2000 GeV, |#eta|>2.196");
	c_norm_C->SaveAs("norm_2000_C_ready.C");
	c_norm_C->SaveAs("norm_2000_C_ready.gif");


	// for saving purpose
	int temp_name;
	TCanvas *c_mod[16]; 	
	for(int i=0; i<16; i++){
		char name[20];
		temp_name = 1+2*i;
		sprintf(name, "c_mod_%i", temp_name);					
		c_mod[i] = new TCanvas(name, name, 1000, 600);
	}

	for(int i=0; i<16; i++){
		c_mod[i]->Divide(2, 1, 0.001, 0.001);
	}

	for(int i=0; i<16; i++){
		temp_mod = 1+2*i;
		//c_mod[i]->cd(1);
		gr_SA[i]->SetTitle("");
		gr_SA[i]->GetXaxis()->SetTitle("#phi [#circ]");
		gr_SA[i]->GetYaxis()->SetTitle("wydajnosc");
		gr_SA[i]->GetHistogram()->SetMinimum(0);
		gr_SA[i]->GetHistogram()->SetMaximum(1);
		//gr_SA[i]->Draw("AP");
		// Tl.DrawLatex(-150, .3, "ATLAS work in progress");
    	// Tl.DrawLatex(-150, .26, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    	// Tl.DrawLatex(-150, .22, "p_{t}>2000 GeV, |#eta|>2.196");
		// c_mod[i]->cd(2);
		// gr_BA[i]->SetTitle("");
		// gr_BA[i]->GetXaxis()->SetTitle("#phi [#circ]");
		// gr_BA[i]->GetYaxis()->SetTitle("wydajnosc");
		// gr_BA[i]->GetHistogram()->SetMinimum(0);
		// gr_BA[i]->GetHistogram()->SetMaximum(1);
		// gr_BA[i]->Draw("AP");
		// Tl.DrawLatex(-150, .9, "ATLAS work in progress");
    	// Tl.DrawLatex(-150, .86, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    	// Tl.DrawLatex(-150, .82, "p_{t}>2000 GeV, |#eta|>2.196");
		char name[20];
		temp_name = 1+2*i;
		sprintf(name, "gr_SA_%i", temp_name);
		gr_SA[i]->SetName(name);
		
	}

	for(int i=0; i<16; i++){
		char name[20];
		temp_name = 1+2*i;
		sprintf(name, "mod_%i_2000_ready.C", temp_name);
		c_mod[i]->SaveAs(name);
		sprintf(name, "mod_%i_2000_ready.gif", temp_name);
		c_mod[i]->SaveAs(name);
	}

	TFile *f1 = new TFile("gr_SA_2000_v2.root", "NEW");
	
	for(int i=0; i<16; i++){
		gr_SA[i]->Write();
	}

	hist_C->Write();
	hist_A->Write();

	h_B[1]->Write();
	h_S[1]->Write();
	h_B[17]->Write();
	h_S[17]->Write();

	f1->Close();

}