void fit_graph_v2(){

	/*
	fitting line to graphs in linear areas (different for each module)
	*/

	TFile *f = new TFile("gr_SA_2000.root", "READ");
	
	// reading from file
	TGraphAsymmErrors *gr[16];
	int temp_mod;

	for(int i=0; i<16; i++){
		temp_mod = 1+2*i;
		char name[20];
		sprintf(name, "gr_SA_%i", temp_mod);
		gr[i] = (TGraphAsymmErrors *) f->Get(name);
		} 

	// module's centre
	float centro[16];
	for(int i=0; i<4; i++){
		centro[i] = 22.5 + 45*i;
		centro[7-i] = (-1)*centro[i];
		centro[i+8] = centro[i];
		centro[15-i] = centro[7-i];
	}
	float phi_min[16];
	float phi_max[16];
	for(int i=0; i<4; i++){
		phi_min[i] = centro[i] - 45;
		phi_max[i] = centro[i] + 45;
	}
	phi_max[3] = -157.5;
	for(int i=0; i<4; i++){
		phi_min[7-i] = -phi_max[i];
		phi_max[7-i] = -phi_min[i];
	}
	for(int i=0; i<4; i++){
		phi_min[i+8] = phi_min[i];
		phi_max[i+8] = phi_max[i];
		phi_min[15-i] = phi_min[7-i];
		phi_max[15-i] = phi_max[7-i];
	}

	// graph with proper points
	TGraphAsymmErrors *gr_fit[16];
	for(int i=0; i<16; i++){
		gr_fit[i] = new TGraphAsymmErrors();
	}

	int n_points = gr[0]->GetN();
	int n_new = 0;

	float diff;
	double x, y;
	double exl, exh, eyl, eyh;

	// rewrite graphs in proper areas
	for(int i=0; i<16; i++){
		n_new = 0;
		for(int n=0; n<n_points; n++){
			gr[i]->GetPoint(n, x, y);
			exl = gr[i]->GetErrorXlow(n);
			exh = gr[i]->GetErrorXhigh(n);
			eyl = gr[i]->GetErrorYlow(n);
			eyh = gr[i]->GetErrorYhigh(n);
			if(i==3 || i==11){
				if(!((x>112.5 && x<180)||(x<-157.5 && x>-180))){
					gr_fit[i]->SetPoint(n_new, x, y);
					gr_fit[i]->SetPointError(n_new, exl, exh, eyl, eyh);
					n_new++;
				}
			}
			else if(i==4 || i==12){
				if(!((x<-112.5 && x>-180) || (x>157.5 && x<180))){
					gr_fit[i]->SetPoint(n_new, x, y);
					gr_fit[i]->SetPointError(n_new, exl, exh, eyl, eyh);
					n_new++;
				}
			}
			else{
				if(!((x>phi_min[i] && x<phi_max[i]))){
					gr_fit[i]->SetPoint(n_new, x, y);
					gr_fit[i]->SetPointError(n_new, exl, exh, eyl, eyh);
					n_new++;
				}
			}
			
		}
	}

    // here: graph is ok, don't touch it

	// calculating mean
	int n_point[16];
	for(int i=0; i<16; i++){
		n_point[i] = gr_fit[i]->GetN();
	}

	double sum[16];
	for(int i=0; i<16; i++){
		sum[i] = 0;
	}

	for(int j=0; j<16; j++)
		{
		for(int i=0; i<n_point[j]; i++){
			gr_fit[j]->GetPoint(i, x, y);
			sum[j]+=y;
		}
	}
	
	for(int i=0; i<16; i++){
		sum[i]/=n_point[i];
	}

	// calculating errors

	double err[16];
	for(int i=0; i<16; i++){
		err[i] = 0;
	}

	float temp_err;

	for(int j=0; j<16; j++){
		for(int i=0; i<n_point[j]; i++){
			temp_err = (gr_fit[j]->GetErrorYhigh(i)+gr_fit[j]->GetErrorYlow(i))/2;
			err[j]+=temp_err*temp_err;
		}
		err[j] = sqrt(err[j]/n_point[j]/n_point[j]);
		
	}

	// mean of fits
	// and error

	float mean = 0;
	float error = 0;

	for(int i=0; i<16; i++){
		mean+=sum[i];
	}
	mean/=16;

	for(int i=0; i<16; i++){
		error+=err[i]*err[i];
	}
	error = sqrt(error/16);

	cout << "pt>2000 MeV" << endl;

	for(int i=0; i<16; i++){
		cout << "modul: " << 2*i+1;
		cout << " fit: " << sum[i] << " blad: " << err[i] << endl;
	}

	cout << "srednio: " << mean << " blad: " << error << endl;

	// for all of them, how it looks

	double mod[16];
	double ex[16];
	for(int i=0; i<16; i++){
		mod[i] = 2*i+1;
		ex[i] = 0;
	}
	double mean_gr[16];
	double err_high[16];
	double err_low[16];
	for(int i=0; i<16; i++){
		mean_gr[i] = mean;
		err_high[i] = mean + error;
		err_low[i] = mean - error; 
	}

	// drawing
	TLatex Tl;
    Tl.SetTextAlign(11);
    Tl.SetTextSize(0.05);
    Tl.SetNDC();

	TCanvas *eff_mod_2000 = new TCanvas("eff_mod_2000", "eff_mod_2000", 1000, 600);
	TGraphErrors *gr_all = new TGraphErrors(16, mod, sum, ex, err);
	gr_all->SetTitle("");
	gr_all->GetXaxis()->SetTitle("numer modulu");
	gr_all->GetYaxis()->SetTitle("wydajnosc");
	gr_all->GetHistogram()->SetMinimum(0);
	gr_all->GetHistogram()->SetMaximum(1);
	gr_all->Draw("AP*");
	TGraph *gr_mean = new TGraph(16, mod, mean_gr);
	gr_mean->SetLineColor(9);
	gr_mean->Draw("L");
	TGraph *gr_err_high = new TGraph(16, mod, err_high);
	gr_err_high->SetLineColor(9);
	gr_err_high->SetLineStyle(2);
	gr_err_high->Draw("L");
	TGraph *gr_err_low = new TGraph(16, mod, err_low);
	gr_err_low->SetLineColor(9);
	gr_err_low->SetLineStyle(2);
	gr_err_low->Draw("L");
	Tl.DrawLatex(0.2, 0.3, "ATLAS work in progress");
    Tl.DrawLatex(.2, .25, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(.2, .20, "p_{t}>2000 GeV, |#eta|>2.196");
	eff_mod_2000->SaveAs("eff_mod_2000.C");
	eff_mod_2000->SaveAs("eff_mod_2000.gif");


}