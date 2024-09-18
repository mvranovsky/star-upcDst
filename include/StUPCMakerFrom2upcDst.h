
#ifndef StUPCMakerFrom2upcDst_h
#define StUPCMakerFrom2upcDst_h

class StUPCMakerFrom2upcDst {

 public:

  StUPCMakerFrom2upcDst(const std::string& in1_name, const std::string& in2_name, const std::string& out_name);
  virtual ~StUPCMakerFrom2upcDst() {}

  Int_t Init() { return 0; }
  Int_t EventLoop();
  Int_t Finish();


 private:
  //void setUpcDstEvent();
  void fillPicoInfo();
  void enableBranchAndSubbranches(TTree* tree, const char* branchName);

  TFile *mInpMain; // main input
  TFile *mInpIdx; // indexed secondary input

  TTree *mTreeMain; // UPC tree from the main input
  TTree *mTreeIdx; // UPC tree from indexed secondary input

  StUPCEvent *mEvtMain; // main input UPC event
  StUPCEvent *mEvtIdx; // indexed secondary input UPC event

  TFile *mOutFile; // output file
  TTree *mTreeOut; // output UPC tree
  StUPCEvent *mEvtOut; // output UPC event

  ClassDef(StUPCMakerFrom2upcDst, 1);

};

#endif

