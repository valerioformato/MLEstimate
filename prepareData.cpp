float tau = 2.2e-6;

void prepareData() {
  gRandom = new TRandom3(0);

  auto cc = new TCanvas();
  cc->Draw();

  auto fun = new TF1("exp", "[0] * exp(-x/[1])", 0, 1e-4);
  fun->SetParameters(1 / tau, tau);
  fun->SetNpx(10000);

  auto hist = new TH1F("hist", ";t (s);Counts", 1000, 0, 1e-4);
  auto tree = new TTree("Triggers", "Triggers");

  float t = 0, dt = 0;
  tree->Branch("time", &t);

  Long64_t nEv = gRandom->Poisson(500);
  cout << "Generating " << nEv << " events" << endl;

  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    t = fun->GetRandom();
    hist->Fill(t);

    tree->Fill();
  }

  cout << "...done" << endl;

  hist->SetStats(0);
  hist->SetMarkerStyle(20);
  hist->GetXaxis()->SetRangeUser(0, 1.1e-5);
  fun->SetLineColor(4);

  hist->Draw("E");
  // fun->SetParameter(0, nEv / tau * hist->GetBinWidth(1));
  // fun->Draw("SAME");

  cc->Print("plots/plot1.pdf");

  auto outfile = new TFile("data/data.root", "recreate");
  outfile->WriteTObject(tree);
  outfile->Close();
}
