#ifndef ITS3_DIGIPARAMS_H
#define ITS3_DIGIPARAMS_H

#include "ITSMFTSimulation/DigiParams.h"

namespace o2 {

namespace its3
{

class DigiParams : public o2::itsmft::DigiParams
{
 public:
  DigiParams() = default;
  ~DigiParams() = default;

  const o2::itsmft::AlpideSimResponse* getAlpSimResponse() const = delete;
  void setAlpSimResponse(const o2::itsmft::AlpideSimResponse* par) = delete;

  const o2::itsmft::AlpideSimResponse* getOBSimResponse() const { return mOBSimResponse; }
  void setOBSimResponse(const o2::itsmft::AlpideSimResponse* response) { mOBSimResponse = response; }

  const o2::itsmft::AlpideSimResponse* getIBSimResponse() const { return mIBSimResponse; }
  void setIBSimResponse(const o2::itsmft::AlpideSimResponse* response) { mIBSimResponse = response; }

  void print() const;

 private:
  const o2::itsmft::AlpideSimResponse* mOBSimResponse = nullptr;  //!< pointer to external response
  const o2::itsmft::AlpideSimResponse* mIBSimResponse = nullptr;  //!< pointer to external response

  ClassDefNV(DigiParams, 1);
};

} // namespace its3

} // namespace o2

#endif