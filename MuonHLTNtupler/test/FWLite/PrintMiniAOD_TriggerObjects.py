# -- reference: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger -- #
# import ROOT in batch mode
import sys
import ROOT
ROOT.gROOT.SetBatch(True)

List_ROOTFile = [
"/u/user/kplee/scratch/ROOTFiles_Test/muonHLTIssues/TripleMu/output1.root",
"/u/user/kplee/scratch/ROOTFiles_Test/muonHLTIssues/TripleMu/output2.root",
"/u/user/kplee/scratch/ROOTFiles_Test/muonHLTIssues/TripleMu/output3.root",
"/u/user/kplee/scratch/ROOTFiles_Test/muonHLTIssues/TripleMu/output4.root",
"/u/user/kplee/scratch/ROOTFiles_Test/muonHLTIssues/TripleMu/output5.root",
]

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerObjects, triggerObjectLabel  = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"
# l1Muons, l1MuonLabel  = Handle("BXVector"), "gmtStage2Digis:Muon"
# l1EGammas, l1EGammaLabel  = Handle("BXVector"), "caloStage2Digis:EGamma"
# l1Jets, l1JetLabel  = Handle("BXVector"), "caloStage2Digis:Jet"
# l1EtSums, l1EtSumLabel  = Handle("BXVector"), "caloStage2Digis:EtSum"
# l1Taus, l1TauLabel  = Handle("BXVector"), "caloStage2Digis:Tau"

events = Events( List_ROOTFile )

for iev,event in enumerate(events):
    event.getByLabel(triggerBitLabel, triggerBits)
    event.getByLabel(triggerObjectLabel, triggerObjects)
    event.getByLabel(triggerPrescaleLabel, triggerPrescales)
    # event.getByLabel(l1MuonLabel, l1Muons)
    # event.getByLabel(l1EGammaLabel, l1EGammas)
    # event.getByLabel(l1JetLabel, l1Jets)
    # event.getByLabel(l1EtSumLabel, l1EtSums)
    # event.getByLabel(l1TauLabel, l1Taus)

    print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    print "\n === TRIGGER PATHS ==="
    names = event.object().triggerNames(triggerBits.product())
    for i in xrange(triggerBits.product().size()):
        print "Trigger ", names.triggerName(i), ", prescale ", triggerPrescales.product().getPrescaleForIndex(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)")

    print "\n === TRIGGER OBJECTS ==="
    for j,to in enumerate(triggerObjects.product()):
        to.unpackPathNames(names);
        print "Trigger object pt %6.2f eta %+5.3f phi %+5.3f  " % (to.pt(),to.eta(),to.phi())
        print "         collection: ", to.collection()
        print "         type ids: ", ", ".join([str(f) for f in to.filterIds()])
        print "         filters: ", ", ".join([str(f) for f in to.filterLabels()])
        pathslast = set(to.pathNames(True))
        print "         paths:   ", ", ".join([("%s*" if f in pathslast else "%s")%f for f in to.pathNames()])

    # print "\n === BARE L1 OBJECTS ==="
    # for where, what in (l1Muons, l1MuonLabel), (l1EGammas, l1EGammaLabel), (l1Jets, l1JetLabel), (l1Taus, l1TauLabel):
    #     shortname = what.split(":")[1]
    #     bxvector = where.product()
    #     for bx in 0,: # xrange(bxvector.getFirstBX(), bxvector.getLastBX()+1): # typically we want only bx=0
    #         for i in xrange(bxvector.size(bx)):
    #             l1obj = bxvector.at(bx,i)
    #             if shortname != "Muon" and l1obj.pt() &lt;= 0.01: continue
    #             print "%-10s  bx %+1d  pt %6.2f eta %+5.3f phi %+5.3f " % (shortname, bx, l1obj.pt(), l1obj.eta(), l1obj.phi())
    # for shortname in ("TotalEt", "TotalHt", "MissingEt", "MissingHt"):
    #     bxvector = l1EtSums.product()
    #     for bx in 0,: # xrange(bxvector.getFirstBX(), bxvector.getLastBX()+1): # typically we want only bx=0
    #         for i in xrange(bxvector.size(bx)):
    #             l1obj = bxvector.at(bx,i)
    #             if l1obj.getType() != getattr(l1obj, "k"+shortname): continue
    #             print "%-10s  bx %+1d  pt %6.2f eta %+5.3f phi %+5.3f " % (shortname, bx, l1obj.pt(), l1obj.eta(), l1obj.phi())
    # if iev > 10: break