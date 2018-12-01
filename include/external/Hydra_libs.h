/*---------------------------------
 * Include hydra classes and
 * algorithms for
 *--------------------------------
 */
#include <hydra/Types.h>
#include <hydra/Vector4R.h>
#include <hydra/PhaseSpace.h>
#include <hydra/Evaluate.h>
#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/FunctionWrapper.h>
#include <hydra/Algorithm.h>
#include <hydra/Tuple.h>
#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/Decays.h>
#include <hydra/functions/WignerDMatrix.h>
#include <hydra/Complex.h>
#include <hydra/functions/BreitWignerLineShape.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/BlattWeisskopfFunctions.h>
#include <hydra/PhaseSpaceIntegrator.h>
#include <hydra/LogLikelihoodFCN.h>
#include <hydra/Parameter.h>
#include <hydra/UserParameters.h>
#include <hydra/Pdf.h>
#include <hydra/AddPdf.h>
#include <hydra/Placeholders.h>
#include <hydra/Random.h>
#include <hydra/Filter.h>
#include <hydra/DenseHistogram.h>
#include <hydra/functions/BreitWignerNR.h>
#include <hydra/functions/Chebychev.h>
#include <hydra/Range.h>
#include <hydra/Distance.h>
#include <hydra/SparseHistogram.h>
#include <hydra/functions/BreitWignerLineShape.h>
#include <hydra/functions/CosHelicityAngle.h>
#include <hydra/functions/ZemachFunctions.h>
