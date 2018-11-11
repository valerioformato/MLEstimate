struct measurement {
  std::string name;
  float value;
  float variance;
  float GetError() { return sqrt(variance); };
  void Print() {
    cout << " - " << name << " = ";
    cout << value << " +- " << GetError();
    cout << "   (variance = " << variance << ")" << endl;
  }
};

measurement fitAnalytical(TTree *tree);
measurement fitScanLL(TTree *tree);
measurement fitMinuit(TTree *tree);
measurement fitMinuit2(TTree *tree);

void fitData() {
  auto infile = std::unique_ptr<TFile>(TFile::Open("data/data.root"));

  auto tree = (TTree *)infile->Get("Triggers");

  auto tauAn = fitAnalytical(tree);
}

measurement fitAnalytical(TTree *tree) {
  float t;
  tree->SetBranchAddress("time", &t);

  auto hist = new TH1F("hist", ";t (s);Counts", 1000, 0, 1e-4);

  measurement tauML;
  tauML.name = "tau";

  Long64_t nEv = tree->GetEntries();
  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    tree->GetEntry(iEv);

    hist->Fill(t);

    tauML.value += t;
  }
  tauML.value /= nEv;
  tauML.variance = (tauML.value * tauML.value) / nEv;

  cout << "Analytical estimate" << endl;
  tauML.Print();

  return tauML;
}
