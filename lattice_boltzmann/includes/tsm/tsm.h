#ifndef __TSM__
#define __TSM__

#include "camera.h"
#include <glm/glm.hpp>
#include <vector>

void computeTsmMatrix(glm::mat4& N_T, const Camera& camera, const glm::mat4& Ml,
                      /*float tsmDistance, (not sure how to use this)*/
                      float percentage, int shadowmapHeight);

void debugTSM(const Camera& camera, const glm::mat4& Ml, float percentage);

#endif //__TSM__
