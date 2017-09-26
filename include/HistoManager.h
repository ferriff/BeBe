#ifndef HistoManager_hh
#define HistoManager_hh

//
// Author: Federico Ferri (Saclay) 2009
//

// HistoManager histos;
// add a template histogram called ZZmass
// histos.addTemplate<TH1F>("ZZmass", new TH1F( "ZZmass", "ZZmass", 100, 50., 150.) );
// [...]
// histos.h<TH1F>("ZZmass", "mass_after_eleID")->Fill(4lMass);>
// [..]
// histos.h<TH1F>("ZZmass", "mass_after_isolation")->Fill(4lMass);
// [...]
// std::string outPlotFile_ = "toto.root";
// histos.save(outPlotFile_.c_str());

#include "TFile.h"

#include <iostream>
#include <map>
#include <string>

#include <cassert>

// Instances of this class should be called prior any TFile opening
// also, possibly for some mischievous ROOT reasons, the save() method
// has to be called to avoid a double-linked corruption (unless there is
// an unspotted bug somewhere here)
class HistoManager
{
        public:
                typedef std::map<std::string, TObject * > Map;

                template <class T> void addTemplate(const char *type, T * templ)
                {
                        assert(m_templates.find(type) == m_templates.end());
                        m_templates.insert(std::make_pair(type, std::move(templ)));
                }

                template <class T> T * h(const char *type, const char *name, T * h = 0)
                {
                        if (h) return h;
                        assert(m_templates.find(type) != m_templates.end());
                        std::string id(type);
                        id.append("_");
                        id.append(name);
                        Map::const_iterator it = m_histos.find(id.c_str());
                        if (it != m_histos.end()) {
                                return (T *)it->second;
                        } else {
                                m_histos[id.c_str()] = (T *)m_templates[type]->Clone(id.c_str());
                                return (T *)m_histos[id.c_str()];
                        }
                }

                void save(const char *fileName = 0)
                {
                        TFile * f(0);
                        TDirectory * d = TDirectory::CurrentDirectory();
                        if (fileName != 0) f = TFile::Open(fileName, "RECREATE");
                        save(f);
                        if (f != 0) f->Close();
                        d->cd();
                        delete f;
                }

                void save(TFile * f)
                {
                        TDirectory * d = TDirectory::CurrentDirectory();
                        assert(f != 0);
                        f->cd();
                        for (auto & it : m_histos) it.second->Write();
                        d->cd();
                }

                void print() const
                {
                        for (auto & it : m_templates) {
                                std::cerr << "Template name: " << it.first
                                          << " object name: " << it.second->GetName()
                                          << " title: " << it.second->GetTitle()
                                          << " pointer: " << &it.second
                                          << "\n";
                        }
                        for (auto & it : m_histos) {
                                std::cerr << "Template name: " << it.first
                                          << " object name: " << it.second->GetName()
                                          << " title: " << it.second->GetTitle()
                                          << " pointer: " << &it.second
                                          << "\n";
                        }
                }

                void clear() { m_histos.clear(); }

        private:
                Map m_templates;
                Map m_histos;
};

#endif
