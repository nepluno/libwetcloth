//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#define GL_SILENCE_DEPRECATION

#include "TwoDSceneRenderer.h"

#include "AttachForce.h"
#include "DistanceFields.h"
#include "MathUtilities.h"
#include "TwoDimensionalDisplayController.h"

const Vector3s vertex_color =
    Vector3s(0.650980392156863, 0.294117647058824, 0.0);
const Vector3s edge_color = Vector3s(0.0, 0.388235294117647, 0.388235294117647);
const Vector3s def_grad_color = Vector3s(0.85, 0.45, 0.0);
const Vector3s sat_color = Vector3s(0.0, 0.0, 1.0);
const Vector3s node_color_x =
    Vector3s(0.650980392156863, 0.1470588235, 0.1470588235);
const Vector3s node_color_y =
    Vector3s(0.1470588235, 0.650980392156863, 0.1470588235);
const Vector3s node_color_z =
    Vector3s(0.1470588235, 0.650980392156863, 0.650980392156863);
const Vector3s node_color_ex = Vector3s(0.451, 0.054, 0.298);
const Vector3s node_color_ey = Vector3s(0.298, 0.451, 0.054);
const Vector3s node_color_ez = Vector3s(0.054, 0.298, 0.451);
const Vector3s node_color_p = Vector3s(0.0, 0.0, 0.0);
const Vector3s node_color_solid_phi = Vector3s(0.85, 0.85, 0.0);
const Vector3s gauss_color =
    Vector3s(0.388235294117647, 0.388235294117647, 0.0);
const Vector3s fluid_color = Vector3s(0, 0, 1);
const Vector3s bucket_color = Vector3s(0.85, 0.85, 0.85);
const Vector3s face_color = Vector3s(0.75, 0.75, 0.75);
const Vector3s face_color2 = Vector3s(0.45, 0.45, 0.45);
const Vector3s attach_color = Vector3s(1.0, 0, 0.);

TwoDSceneRenderer::TwoDSceneRenderer(const TwoDScene& scene) {
  m_info.render_particles = true;
  m_info.render_vertices = false;
  m_info.render_gauss = false;
  m_info.render_particle_velocity = false;
  m_info.render_vertice_velocity = false;
  m_info.render_gauss_velocity = false;
  m_info.render_cloth = true;
  m_info.render_yarn = true;
  m_info.render_levelset = true;
  m_info.render_spring = true;
  m_info.render_cohesion = false;
  m_info.render_deformation_gradient_length = 0.0;
  m_info.render_buckets = false;
  m_info.render_velocity_length = 10.0;

  m_info.render_nodes = RenderInfo::RI_NV_NONE;
  m_info.render_face_centers = RenderInfo::RI_FCV_NONE;
  m_info.render_edge_centers = RenderInfo::RI_ECV_NONE;
  m_info.render_cell_centers = RenderInfo::RI_CCV_NONE;

  const auto& groups = scene.getParticleGroup();
  if (groups.size() > 0) {
    const int num_group = *std::max_element(groups.begin(), groups.end()) + 1;
    m_group_colors.resize(num_group);

    for (int i = 0; i < num_group; ++i) {
      m_group_colors[i] = Vector3s::Constant((scalar)i / (scalar)num_group);
    }
  }
}

void TwoDSceneRenderer::updateParticleSimulationState(const TwoDScene& scene) {}

RenderInfo& TwoDSceneRenderer::getRenderInfo() { return m_info; }

const RenderInfo& TwoDSceneRenderer::getRenderInfo() const { return m_info; }

void TwoDSceneRenderer::renderParticleSimulation(const TwoDScene& scene,
                                                 const scalar& dt) {
#ifdef RENDER_ENABLED
  const VectorXs& x = scene.getX();
  const VectorXs& rest_x = scene.getRestPos();
  const VectorXs& gx = scene.getGaussX();

  const VectorXs& v = scene.getV();
  const VectorXs& gv = scene.getGaussV();

  const VectorXs& vol = scene.getVol();
  const VectorXs& fvol = scene.getFluidVol();
  const VectorXs& fv = scene.getFluidV();

  const MatrixXi& faces = scene.getFaces();
  const MatrixXs& fe = scene.getGaussFe();
  const int num_faces = faces.rows();

  const std::vector<std::vector<RayTriInfo> >& intersections =
      scene.getIntersections();
  const int num_gauss = scene.getNumGausses();

  const scalar dx = scene.getCellSize();

  assert(x.size() % 4 == 0);
  assert(4 * scene.getNumParticles() == x.size());

  // Render faces
  if (m_info.render_cloth) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < num_faces; ++i) {
      const auto& f = faces.row(i);
      for (int r = 0; r < 3; ++r) {
        const scalar sat =
            mathutils::clamp(fvol(f(r)) / std::max(1e-12, vol(f(r))), 0.0, 1.0);
        const Vector3s c = face_color * (1.0 - sat) + fluid_color * sat;
        glColor3d(c(0), c(1), c(2));
        glVertex3dv(x.segment<3>(f(r) * 4).data());
      }
    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < num_faces; ++i) {
      const auto& f = faces.row(i);
      for (int r = 0; r < 3; ++r) {
        const scalar sat =
            mathutils::clamp(fvol(f(r)) / std::max(1e-12, vol(f(r))), 0.0, 1.0);
        const Vector3s c =
            (face_color2 * (1.0 - sat) + fluid_color * sat) * 0.85 +
            face_color2 * 0.15;
        glColor3d(c(0), c(1), c(2));
        glVertex3dv(x.segment<3>(f(r) * 4).data());
      }
    }
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  if (m_info.render_yarn) {
    // Render edges
    const MatrixXi& edges = scene.getEdges();
    const auto& groups = scene.getParticleGroup();

    glLineWidth(3.0);
    glEnable(GL_DEPTH_TEST);
    glBegin(GL_LINES);

    for (int i = 0; i < edges.rows(); ++i) {
      const scalar sat0 = mathutils::clamp(
          fvol(edges(i, 0)) / std::max(1e-12, vol(edges(i, 0))), 0.0, 1.0);
      const scalar sat1 = mathutils::clamp(
          fvol(edges(i, 1)) / std::max(1e-12, vol(edges(i, 1))), 0.0, 1.0);

      const Vector3s c0 =
          m_group_colors[groups[edges(i, 0)] % m_group_colors.size()] *
              (1.0 - sat0) +
          fluid_color * sat0;
      const Vector3s c1 =
          m_group_colors[groups[edges(i, 1)] % m_group_colors.size()] *
              (1.0 - sat1) +
          fluid_color * sat1;

      glColor3dv(c0.data());
      glVertex3dv(x.segment<3>(4 * edges(i, 0)).data());

      glColor3dv(c1.data());
      glVertex3dv(x.segment<3>(4 * edges(i, 1)).data());

      //        renderSweptEdge( x.segment<3>(4 * edges(i, 0)), x.segment<3>(4 *
      //        edges(i, 1)), (scene.getRadius()[edges(i, 0)] +
      //        scene.getRadius()[edges(i, 1)]) * 0.5, c0, c1 );
    }
    glDisable(GL_DEPTH_TEST);
    glEnd();
  }

  if (m_info.render_spring) {
    // Render Spring
    const std::vector<std::shared_ptr<AttachForce> >& attaches =
        scene.getAttachForces();

    glColor3dv(attach_color.data());
    glBegin(GL_LINES);
    for (const std::shared_ptr<AttachForce>& force : attaches) {
      if (force->getKs() == 0.0) continue;

      glVertex3dv(x.segment<3>(force->getParticleIndex() * 4).data());
      glVertex3dv(rest_x.segment<3>(force->getParticleIndex() * 4).data());
    }
    glEnd();
  }

  if (m_info.render_gauss) {
    // render Gauss
    glColor3dv(gauss_color.data());
    glPointSize(5);
    glBegin(GL_POINTS);
    for (int i = 0; i < scene.getNumGausses(); ++i) {
      glVertex3dv(gx.segment<3>(4 * i).data());
    }
    glEnd();
  }

  if (m_info.render_buckets) {
    const Sorter& bucket = scene.getParticleBuckets();
    const Vector3s& min_corner = scene.getBucketMinCorner();
    const scalar bl = scene.getBucketLength();

    glColor4d(0.5, 0.5, 0.5, 0.25);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINES);

    for (int k = 0; k <= bucket.nk; ++k)
      for (int j = 0; j <= bucket.nj; ++j) {
        Vector3s s = Vector3s(0.0, j, k) * bl + min_corner;
        Vector3s e = Vector3s(bucket.ni, j, k) * bl + min_corner;

        glVertex3dv(s.data());
        glVertex3dv(e.data());
      }

    for (int k = 0; k <= bucket.nk; ++k)
      for (int i = 0; i <= bucket.ni; ++i) {
        Vector3s s = Vector3s(i, 0.0, k) * bl + min_corner;
        Vector3s e = Vector3s(i, bucket.nj, k) * bl + min_corner;

        glVertex3dv(s.data());
        glVertex3dv(e.data());
      }

    for (int j = 0; j <= bucket.nj; ++j)
      for (int i = 0; i <= bucket.ni; ++i) {
        Vector3s s = Vector3s(i, j, 0.0) * bl + min_corner;
        Vector3s e = Vector3s(i, j, bucket.nk) * bl + min_corner;

        glVertex3dv(s.data());
        glVertex3dv(e.data());
      }

    glEnd();

    glDisable(GL_BLEND);
  }

  if (m_info.render_cohesion) {
    glBegin(GL_LINES);
    for (int i = 0; i < num_gauss; ++i) {
      const std::vector<RayTriInfo>& inters = intersections[i];
      for (auto& info : inters) {
        const Vector3s c =
            renderingutils::interpolateColor(info.volume_frac, 0.0, 1.0);
        glColor3dv(c.data());

        glVertex3dv(gx.segment<3>(i * 4).data());
        glVertex3dv(info.end.data());
      }
    }
    glEnd();
  }

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBegin(GL_POINTS);
  switch (m_info.render_nodes) {
    case RenderInfo::RI_NV_CONSTANT: {
      glColor4d(node_color_solid_phi(0), node_color_solid_phi(1),
                node_color_solid_phi(2), 0.8);
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const int num_node_sphi = scene.getNumNodes(i);
        for (int j = 0; j < num_node_sphi; ++j) {
          glVertex3dv(scene.getNodePosSolidPhi(i, j).data());
        }
      }
      break;
    }
    case RenderInfo::RI_NV_SOLID_PHI: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_sphi = scene.getNodeSolidPhi()[i];
        const int num_node_sphi = scene.getNumNodes(i);
        for (int j = 0; j < num_node_sphi; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_sphi[j], -3.0 * dx,
                                                        3.0 * dx);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosSolidPhi(i, j).data());
        }
      }
      break;
    }
    default:
      break;
  }

  switch (m_info.render_face_centers) {
    case RenderInfo::RI_FCV_CONSTANT: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_x(0), node_color_x(1), node_color_x(2), 0.8);
        const int num_node_x = scene.getNumNodes(i);
        for (int j = 0; j < num_node_x; ++j) {
          glVertex3dv(scene.getNodePosX(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_y(0), node_color_y(1), node_color_y(2), 0.8);
        const int num_node_y = scene.getNumNodes(i);
        for (int j = 0; j < num_node_y; ++j) {
          glVertex3dv(scene.getNodePosY(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_z(0), node_color_z(1), node_color_z(2), 0.8);
        const int num_node_z = scene.getNumNodes(i);
        for (int j = 0; j < num_node_z; ++j) {
          glVertex3dv(scene.getNodePosZ(i, j).data());
        }
      }
      break;
    }
    case RenderInfo::RI_FCV_SOLID_VOL: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_solid_psi = scene.getNodePsiX()[i];
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_solid_psi[j]);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosX(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_solid_psi = scene.getNodePsiY()[i];
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_solid_psi[j]);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosY(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_solid_psi = scene.getNodePsiZ()[i];
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_solid_psi[j]);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosZ(i, j).data());
        }
      }
      break;
    }
    case RenderInfo::RI_FCV_LIQUID_VOL: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_sat = scene.getNodeSaturationX()[i];
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_sat[j]);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosX(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_sat = scene.getNodeSaturationY()[i];
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_sat[j]);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosY(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        const VectorXs& node_sat = scene.getNodeSaturationZ()[i];
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_sat[j]);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosZ(i, j).data());
        }
      }
      break;
    }
    default:
      break;
  }

  switch (m_info.render_edge_centers) {
    case RenderInfo::RI_ECV_CONSTANT: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_ex(0), node_color_ex(1), node_color_ex(2), 0.8);
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          glVertex3dv(scene.getNodePosX(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_ey(0), node_color_ey(1), node_color_ey(2), 0.8);
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          glVertex3dv(scene.getNodePosY(i, j).data());
        }
      }

      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_ez(0), node_color_ez(1), node_color_ez(2), 0.8);
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          glVertex3dv(scene.getNodePosZ(i, j).data());
        }
      }
      break;
    }
    default:
      break;
  }

  switch (m_info.render_cell_centers) {
    case RenderInfo::RI_CCV_CONSTANT: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_p(0), node_color_p(1), node_color_p(2), 0.8);
        const int num_nodes = scene.getNumNodes(i);
        for (int j = 0; j < num_nodes; ++j) {
          glVertex3dv(scene.getNodePosP(i, j).data());
        }
      }

      break;
    }
    case RenderInfo::RI_CCV_LIQUID_PHI: {
      for (int i = 0; i < scene.getNumBuckets(); ++i) {
        glColor4d(node_color_p(0), node_color_p(1), node_color_p(2), 0.8);
        const int num_nodes = scene.getNumNodes(i);
        const VectorXs& node_phi = scene.getNodeLiquidPhi()[i];
        for (int j = 0; j < num_nodes; ++j) {
          Vector3s c = renderingutils::interpolateColor(node_phi[j], -3.0 * dx,
                                                        3.0 * dx);
          glColor4d(c(0), c(1), c(2), 0.8);
          glVertex3dv(scene.getNodePosP(i, j).data());
        }
      }

      break;
    }
    default:
      break;
  }
  glEnd();
  glDisable(GL_BLEND);

  // render fluid particles
  if (m_info.render_particles) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    const std::vector<int>& fluid_indices = scene.getFluidIndices();
    glPointSize(4.0);
    glColor4d(0.0, 0.0, 1.0, 0.1);
    glBegin(GL_POINTS);
    for (int i : fluid_indices) {
      const Vector3s& px = x.segment<3>(i * 4);
      glVertex3dv(px.data());
    }
    glEnd();
    glDisable(GL_BLEND);
  }

  // render vertices
  if (m_info.render_vertices) {
    const int num_soft_elasto = scene.getNumSoftElastoParticles();
    glPointSize(4.0);
    glColor3dv(vertex_color.data());
    glBegin(GL_POINTS);
    for (int i = 0; i < num_soft_elasto; ++i) {
      const Vector3s& px = x.segment<3>(i * 4);
      glVertex3dv(px.data());
    }
    glEnd();
  }

  if (m_info.render_vertice_velocity) {
    glColor4d(vertex_color(0), vertex_color(1), vertex_color(2), 0.25);
    const int num_soft_elasto = scene.getNumSoftElastoParticles();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINES);
    for (int i = 0; i < num_soft_elasto; ++i) {
      const Vector3s& px = x.segment<3>(i * 4);
      glVertex3dv(px.data());
      const Vector3s pxn =
          x.segment<3>(i * 4) +
          v.segment<3>(i * 4) * dt * m_info.render_velocity_length;
      glVertex3dv(pxn.data());
    }
    glEnd();
    glDisable(GL_BLEND);
  }

  if (m_info.render_gauss_velocity) {
    glColor4d(gauss_color(0), gauss_color(1), gauss_color(2), 0.25);
    const int num_gauss = scene.getNumGausses();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINES);
    for (int i = 0; i < num_gauss; ++i) {
      const Vector3s& px = gx.segment<3>(i * 4);
      glVertex3dv(px.data());
      const Vector3s pxn =
          gx.segment<3>(i * 4) +
          gv.segment<3>(i * 4) * dt * m_info.render_velocity_length;
      glVertex3dv(pxn.data());
    }
    glEnd();
    glDisable(GL_BLEND);
  }

  if (m_info.render_particle_velocity) {
    glColor4d(fluid_color(0), fluid_color(1), fluid_color(2), 0.25);
    const std::vector<int>& fluid_indices = scene.getFluidIndices();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINES);
    for (int i : fluid_indices) {
      const Vector3s& px = x.segment<3>(i * 4);
      glVertex3dv(px.data());
      const Vector3s pxn =
          x.segment<3>(i * 4) +
          fv.segment<3>(i * 4) * dt * m_info.render_velocity_length;
      glVertex3dv(pxn.data());
    }
    glEnd();
    glDisable(GL_BLEND);
  }

  if (m_info.render_deformation_gradient_length > 0.0) {
    glColor4d(def_grad_color(0), def_grad_color(1), def_grad_color(2), 0.25);
    const int num_gauss = scene.getNumGausses();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINES);
    for (int i = 0; i < num_gauss; ++i) {
      const Vector3s& px = gx.segment<3>(i * 4);
      glVertex3dv(px.data());
      const Vector3s px0 = px + fe.block<3, 1>(i * 3, 0) *
                                    m_info.render_deformation_gradient_length;
      glVertex3dv(px0.data());
      glVertex3dv(px.data());
      const Vector3s px1 = px + fe.block<3, 1>(i * 3, 1) *
                                    m_info.render_deformation_gradient_length;
      glVertex3dv(px1.data());
      glVertex3dv(px.data());
      const Vector3s px2 = px + fe.block<3, 1>(i * 3, 2) *
                                    m_info.render_deformation_gradient_length;
      glVertex3dv(px2.data());
    }
    glEnd();
    glDisable(GL_BLEND);
  }

  // render mesh
  if (m_info.render_levelset) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POLYGON_STIPPLE);
    const std::vector<std::shared_ptr<DistanceField> >& fields =
        scene.getGroupDistanceField();
    for (auto& ptr : fields) {
      if (ptr->usage == DFU_SOLID) {
        glColor4d(0.0, 0.0, 0.0, 0.02);
      } else if (ptr->usage == DFU_TERMINATOR) {
        glColor4d(0.0, 0.0, 1.0, 0.02);
      } else {
        continue;
      }

      for (int i = 0; i < 2; ++i) {
        if (i == 1) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        ptr->render([](const std::vector<Vector3s>& vertices,
                       const std::vector<Vector3i>& indices,
                       const Eigen::Quaternion<scalar>& rot,
                       const Vector3s& center, const scalar&) {
          glPushMatrix();
          glTranslated(center(0), center(1), center(2));
          Eigen::AngleAxis<scalar> rotaa(rot);
          glRotated(rotaa.angle() * 180.0 / M_PI, rotaa.axis()(0),
                    rotaa.axis()(1), rotaa.axis()(2));

          glBegin(GL_TRIANGLES);
          for (const Vector3i& idx : indices) {
            for (int r = 0; r < 3; ++r) {
              glVertex3dv(vertices[idx(r)].data());
            }
          }

          glEnd();
          glPopMatrix();
        });
      }

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    glDisable(GL_BLEND);
    glDisable(GL_POLYGON_STIPPLE);
  }
#endif
}
