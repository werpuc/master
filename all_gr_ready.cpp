void all_gr_ready(){

    /*
    macro for analysing data from both discs
    comparison between automating fitting and its errors and those calculated by me
    */

	TFile *f = new TFile("gr_SA_2000_all.root", "READ");
	
	TH1F *all_S_A = (TH1F *) f->Get("all_S_A");
	TH1F *all_norm_A = (TH1F *) f->Get("all_norm_A");

	TH1F *all_S_C = (TH1F *) f->Get("all_S_C");
	TH1F *all_norm_C = (TH1F *) f->Get("all_norm_C");

	TGraphAsymmErrors *gr_A = new TGraphAsymmErrors();
	gr_A->BayesDivide(all_S_A, all_norm_A);
	gr_A->GetHistogram()->SetMinimum(0);
	gr_A->GetHistogram()->SetMaximum(1.02);
	gr_A->SetTitle("");
	gr_A->GetXaxis()->SetTitle("#phi [#circ]");
	gr_A->GetYaxis()->SetTitle("wydajnosc");

	TGraphAsymmErrors *gr_C = new TGraphAsymmErrors();
	gr_C->BayesDivide(all_S_C, all_norm_C);
	gr_C->GetHistogram()->SetMinimum(0);
	gr_C->GetHistogram()->SetMaximum(1.02);
	gr_C->SetTitle("");
	gr_C->GetXaxis()->SetTitle("#phi [#circ]");
	gr_C->GetYaxis()->SetTitle("wydajnosc");

	TLatex Tl;
    Tl.SetTextAlign(11);
    Tl.SetTextSize(0.05);
    Tl.SetNDC();

    // fitting
    TGraphAsymmErrors *gr_fit_A = new TGraphAsymmErrors();
    TGraphAsymmErrors *gr_fit_C = new TGraphAsymmErrors();
    int n_new = 0;

    double x, y;
    double exla, exha, eyla, eyha;
    double exlc, exhc, eylc, eyhc;

    for(int i=0; i<90; i++){
    	gr_A->GetPoint(i, x, y);
    	exla = gr_A->GetErrorXlow(i);
        exha = gr_A->GetErrorXhigh(i);
        eyla = gr_A->GetErrorYlow(i);
        eyha = gr_A->GetErrorYhigh(i);
    	if(x<-45 || x>45){
    		gr_fit_A->SetPoint(n_new, x, y);
    		gr_fit_A->SetPointError(n_new, exla, exha, eyla, eyha);
    		n_new++;
    	}
    }
    n_new = 0;
    for(int i=0; i<90; i++){
    	gr_C->GetPoint(i, x, y);
    	exlc = gr_C->GetErrorXlow(i);
        exhc = gr_C->GetErrorXhigh(i);
        eylc = gr_C->GetErrorYlow(i);
        eyhc = gr_C->GetErrorYhigh(i);
    	if(x<-45 || x>45){
    		gr_fit_C->SetPoint(n_new, x, y);
    		gr_fit_C->SetPointError(n_new, exlc, exhc, eylc, eyhc);
    		n_new++;
    	}
    }
    TF1 *f1_A = new TF1("f1_A", "[0]", -180, 180);
    TF1 *f1_C = new TF1("f1_C", "[0]", -180, 180);
    f1_A->SetLineStyle(2);
    f1_C->SetLineStyle(2);
    gr_fit_A->Fit("f1_A", "R");
    gr_fit_C->Fit("f1_C", "R");

// manual fitting
    int n_point_A = gr_fit_A->GetN();
    int n_point_C = gr_fit_C->GetN();
    float temp_err;
    float err_fit_a, err_fit_c;
    err_fit_a = err_fit_c = 0;

    for(int i=0; i<n_point_A; i++){
            temp_err = (gr_fit_A->GetErrorYhigh(i)+gr_fit_A->GetErrorYlow(i))/2;
            err_fit_a+=temp_err*temp_err;
        }
    err_fit_a = sqrt(err_fit_a/n_point_A/n_point_A);
    
    for(int i=0; i<n_point_C; i++){
            temp_err = (gr_fit_C->GetErrorYhigh(i)+gr_fit_C->GetErrorYlow(i))/2;
            err_fit_c+=temp_err*temp_err;
        }
    err_fit_c = sqrt(err_fit_c/n_point_C/n_point_C);

    cout << "fit reczny: " << endl;
    cout << "blad A: " << err_fit_a << endl;
    cout << "blad C: " << err_fit_c << endl;



// drawing
    TGraph *line_a_A = new TGraph();
    TGraph *line_b_A = new TGraph();
    TGraph *line_a_C = new TGraph();
    TGraph *line_b_C = new TGraph();

    int n_new_a, n_new_b;
    n_new_a = n_new_b = 0;
    for(int i=0; i<gr_fit_A->GetN(); i++){
    	gr_fit_A->GetPoint(i, x, y);
    	if(x< 0){
    		line_a_A->SetPoint(n_new_a, x, f1_A->GetParameter(0));
    		n_new_a++;
    	}
    	else{
    		line_b_A->SetPoint(n_new_b, x, f1_A->GetParameter(0));
    		n_new_b++;
    	}
    }
    n_new_a = n_new_b = 0;
    for(int i=0; i<gr_fit_C->GetN(); i++){
    	gr_fit_C->GetPoint(i, x, y);
    	if(x< 0){
    		line_a_C->SetPoint(n_new_a, x, f1_C->GetParameter(0));
    		n_new_a++;
    	}
    	else{
    		line_b_C->SetPoint(n_new_b, x, f1_C->GetParameter(0));
    		n_new_b++;
    	}
    }

    line_a_A->SetLineWidth(2);
    line_a_A->SetLineColor(2);
    line_b_A->SetLineWidth(2);
    line_b_A->SetLineColor(2);
    line_a_C->SetLineWidth(2);
    line_a_C->SetLineColor(2);
    line_b_C->SetLineWidth(2);
    line_b_C->SetLineColor(2);

// end of fitting
	/*TCanvas *c_all_disc = new TCanvas("c_all_disc", "c_all_disc", 1000, 600);
	c_all_disc->Divide(1, 2);
	c_all_disc->cd(1);
	gr_A->Draw("AP");
	gr_fit_A->Draw("P");
	line_a_A->Draw("L");
    line_b_A->Draw("L");
	Tl.DrawLatex(0.2, 0.3, "ATLAS work in progress");
    Tl.DrawLatex(.2, .25, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(.2, .20, "p_{t}>2000 GeV, |#eta|>2.196");  
	c_all_disc->cd(2);
	gr_C->Draw("AP");
	gr_fit_C->Draw("P");
	line_a_C->Draw("L");
    line_b_C->Draw("L");
	Tl.DrawLatex(0.2, 0.3, "ATLAS work in progress");
    Tl.DrawLatex(.2, .25, "p+Pb, Pb+p, #sqrt{s} = 5.02TeV");
    Tl.DrawLatex(.2, .20, "p_{t}>2000 GeV, |#eta|>2.196");  
    c_all_disc->SaveAs("eff_all_disc_2000.C");
    c_all_disc->SaveAs("eff_all_disc_2000.gif");
*/


// once again, different method
    float sum_a = 0;
    float sum_a_norm = 0;
    float sum_a_e = 0;
    float sum_a_norm_e = 0;
    float sum_c = 0;
    float sum_c_norm = 0;
    float sum_c_e = 0;
    float sum_c_norm_e = 0;

    for(int i=42; i<50; i++){
    	sum_a+=all_S_A->GetBinContent(i);
    	sum_a_norm+=all_norm_A->GetBinContent(i);
    	//sum_a_e+=pow(all_S_A->GetBinError(i), 2);
    	//sum_a_norm_e+=pow(all_norm_A->GetBinError(i), 2);
    	sum_c+=all_S_C->GetBinContent(i);
    	sum_c_norm+=all_norm_C->GetBinContent(i);
    	//sum_c_e+=pow(all_S_C->GetBinError(i), 2);
    	//sum_c_norm_e+=pow(all_norm_C->GetBinError(i), 2);
    }


    sum_a_e = sqrt(sum_a);
    sum_a_norm_e = sqrt(sum_a_norm);
    sum_c_e = sqrt(sum_c);
    sum_c_norm_e = sqrt(sum_c_norm);

    double sum_err_a, sum_err_c;

    sum_err_a = sqrt(pow( sum_a_e/sum_a_norm,2)+pow(sum_a_norm_e*(sum_a/sum_a_norm/sum_a_norm),2));
    sum_err_c = sqrt(pow( sum_c_e/sum_c_norm,2)+pow(sum_c_norm_e*(sum_c/sum_c_norm/sum_c_norm),2));

    cout << "a: " << sum_a/sum_a_norm << endl;
    cout << "blad: " << sum_err_a << endl;
    cout << "c: " << sum_c/sum_c_norm << endl;
    cout << "blad: " << sum_err_c << endl;

    cout << "************" << endl;
    cout << sum_c << endl;
    cout << sum_c_norm << endl;

    TH1F *h_n_a = new TH1F("h_n_a", "h_n_a", 1, 0, 1);
    TH1F *h_n_syg_a = new TH1F("h_n_syg_a", "h_n_syg_a", 1, 0, 1);
    TH1F *h_n_c = new TH1F("h_n_c", "h_n_c", 1, 0, 1);
    TH1F *h_n_syg_c = new TH1F("h_n_syg_c", "h_n_syg_c", 1, 0, 1);

    h_n_a->SetBinContent(1, sum_a_norm);
    h_n_syg_a->SetBinContent(1, sum_a);
    h_n_c->SetBinContent(1, sum_c_norm);
    h_n_syg_c->SetBinContent(1, sum_c);

    TGraphAsymmErrors *gr_a_one = new TGraphAsymmErrors();
    TGraphAsymmErrors *gr_c_one = new TGraphAsymmErrors();

    gr_a_one->BayesDivide(h_n_syg_a, h_n_a);
    gr_c_one->BayesDivide(h_n_syg_c, h_n_c);    

    cout << "a" << endl;
    gr_a_one->GetPoint(0, x, y);
    cout << y << endl;
    cout << gr_a_one->GetErrorYhigh(0) << endl;
    cout << gr_a_one->GetErrorYlow(0) << endl;

    cout << "c" << endl;
    gr_c_one->GetPoint(0, x, y);
    cout << y << endl;
    cout << gr_c_one->GetErrorYhigh(0) << endl;
    cout << gr_c_one->GetErrorYlow(0) << endl;




}
