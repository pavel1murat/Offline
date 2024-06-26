
#include <iostream>

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetector.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetectorMaker.hh"

using namespace std;

namespace mu2e {


  Mu2eDetector::ptr_t Mu2eDetectorMaker::fromFcl(
      Mu2eMaterial::cptr_t mat_p, Tracker::cptr_t trk_p) {

    Mu2eDetector::ptr_t ptr = make_shared<Mu2eDetector>();

    // loop over Planes
    Mu2eMaterial const& material = *mat_p;
    Tracker const& tracker = *trk_p;
    for ( size_t i=0; i!= tracker.nPlanes(); ++i){
      StrawId sid(i,0,0);
      if(tracker.planeExists(sid)){
        const auto& plane = tracker.getPlane(i);
        // loop over panels
        for(auto panel_p : plane.getPanels()){
          auto& panel = *panel_p;
          // loop over straws
          for (const auto& straw : panel.getStrawPointers()) {
            // build the straw elements from this
            // have to strip const because thing inside BTrk are non-const
            auto temp = const_cast<DetStrawType*>(material.strawType());
            DetStrawElem* elem = new DetStrawElem(temp,tracker,straw->id());
            // push this into the map
            ptr->_strawmap[straw->id()] = elem;
          } // straws
        } // panels
      } // if exists
    } // planes

    return ptr;
  }


} // namespace mu2e
