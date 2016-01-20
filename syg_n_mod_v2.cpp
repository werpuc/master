void syg_n_mod_v2(){

	/*
	how many modules gives signal id only one track is found
	pt>2000
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

	int j; 
	
	int tmp; // dla odwolan
	float temp; // dla odwolan

	int trkA[2]; // numery sladow na dysk w przedziale
	int trkC[2];
	int nA, nC;

	int trk_a, trk_c; 

	float eta, phi; 
	float eta_gr; 

	float pt_cut = 2000;

	int mod_a, mod_c;

	TH1F *h_n_mod_500 = new TH1F("h_n_mod_500", " ; n_{mod}; p-two n_{mod} dajacych sygnal", 9, -0.5, 8.5);
	TH1F *h_n_mod_2000 = new TH1F("h_n_mod_2000", " ; n_{mod}; p-two n_{mod} dajacych sygnal", 9, -0.5, 8.5);

	cout << "obliczenia: start" << endl;

	for(j=0; ; j++){
		int st = hichain.GetEntry(j);
		if(st<1) break;	
		// ************
		//odwolania
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
			if(eta>2){
				if(nA==0){
					trkA[0] = i;
					nA++;
				}
				else {
					trkA[1] = i;
					nA++;
				}
			}
			if(eta<-2){
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
		if(trkA[0]==-1 && trkC[0]==-1) continue;
		if(trkA[1]!=-1){
			trkA[0] = -1;
			trkA[1] = -1;
		}
		if(trkC[1]!=-1){
			trkC[0] = -1;
			trkC[1] = -1;
		}
		if(trkA[0]==-1 && trkC[0]==-1) continue;
		// mamy: wybrane przypadki z tylko jednym śladem na dysk
		trk_a = trk_c = -1;
		trk_a = trkA[0];
		trk_c = trkC[0];
		if(j%1000==0) cout << "j: " << j << endl;
		mod_a = mod_c = 0;		
		if(trkA[0]!=-1){
			if(hichain.trk_pt->at(trkA[0])>2000){
				for(int i=16; i<32; i++){
					if(i%2==0) continue;
					if(hichain.mb_E->at(i)>0.1){
						mod_a++;
					}
				}
				h_n_mod_2000->Fill(mod_a);
			}			
		}
		if(trkC[0]!=-1){
			if(hichain.trk_pt->at(trkC[0])>2000){
				for(int i=0; i<16; i++){
					if(i%2==0) continue;
					if(hichain.mb_E->at(i)>0.1){
						mod_c++;
					}
				}
				h_n_mod_2000->Fill(mod_c);
			}			
		}
		mod_a = mod_c = 0;		
		if(trkA[0]!=-1){
			if(hichain.trk_pt->at(trkA[0])>500){
				for(int i=17; i<32; i=i+2){
					if(hichain.mb_E->at(i)>0.1){
						mod_a++;
					}
				}
				h_n_mod_500->Fill(mod_a);
			}			
		}
		if(trkC[0]!=-1){
			if(hichain.trk_pt->at(trkC[0])>500){
				for(int i=1; i<16; i=i+2){
					if(hichain.mb_E->at(i)>0.1){
						mod_c++;
					}
				}
				h_n_mod_500->Fill(mod_c);
			}			
		}			
	}

	cout << "obliczenia: stop " << endl;

	TLatex Tl;
    Tl.SetTextAlign(11);
    Tl.SetTextSize(0.05);
    Tl.SetNDC();

    double n_scal_500, n_scal_2000;
    n_scal_500 = h_n_mod_500->GetEntries();
    n_scal_2000 = h_n_mod_2000->GetEntries();

    cout << "scal 500: " << n_scal_500 << endl;
    cout << "scal 2000: " << n_scal_2000 << endl;

    h_n_mod_500->Scale(1/n_scal_500);
    h_n_mod_2000->Scale(1/n_scal_2000);

    h_n_mod_2000->SetMinimum(0);
    h_n_mod_2000->SetMaximum(0.2);
    h_n_mod_500->SetMinimum(0);
    h_n_mod_500->SetMaximum(0.2);

	TCanvas *c_n_mod = new TCanvas("c_n_mod", "c_n_mod", 1000, 600);
	c_n_mod->Divide(1, 2);
	c_n_mod->cd(1);
	h_n_mod_500->SetLineWidth(2);
	h_n_mod_500->GetXaxis()->SetLabelSize(0.05);
	h_n_mod_500->GetYaxis()->SetLabelSize(0.05);
	h_n_mod_500->SetTitleSize(0.05, "X");
	h_n_mod_500->SetTitleSize(0.05, "Y");
    h_n_mod_500->Draw();
    Tl.DrawLatex(0.6, 0.3, "ATLAS work in progress");
    Tl.DrawLatex(.6, .25, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(.6, .2, "p_{t}>500 GeV");  
    c_n_mod->cd(2);
    h_n_mod_2000->SetLineWidth(2);
    h_n_mod_2000->GetXaxis()->SetLabelSize(0.05);
	h_n_mod_2000->GetYaxis()->SetLabelSize(0.05);
    h_n_mod_2000->SetTitleSize(0.05, "X");
	h_n_mod_2000->SetTitleSize(0.05, "Y");
	h_n_mod_2000->Draw();
	Tl.DrawLatex(0.6, 0.3, "ATLAS work in progress");
    Tl.DrawLatex(0.6, 0.25, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(0.6, 0.2, "p_{t}>2000 GeV");  
    
	c_n_mod->SaveAs("syg_n_mod_v2.png");

	
}