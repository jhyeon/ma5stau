#ifndef analysis_ATLAS_SUSY_2018_04_h
#define analysis_ATLAS_SUSY_2018_04_h

#include "SampleAnalyzer/Interfaces/root/RootMainHeaders.h"
#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class ATLAS_SUSY_2018_04 : public AnalyzerBase
{
  INIT_ANALYSIS(ATLAS_SUSY_2018_04,"ATLAS_SUSY_2018_04")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
};
}

#endif
