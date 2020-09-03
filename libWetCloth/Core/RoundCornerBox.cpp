//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "RoundCornerBox.h"

RoundCornerBox::RoundCornerBox(int N, const Vector3s& b, const double& radius) {
  int N_edge = 2 * (N + 1);
  m_N_edge = N_edge;
  m_index_to_verts.resize(N_edge * N_edge * N_edge, -1);
  m_radius = radius;
  double dx = radius / (double)N;

  const double sign[2] = {-1.0, 1.0};
  const int ks[2] = {0, N * 2 + 1};

  // xy-planes
  for (int kidx = 0; kidx < 2; ++kidx) {
    const int k = ks[kidx];
    Vector3s origin =
        Vector3s(-b(0) - radius, -b(1) - radius, (b(2) + radius) * sign[kidx]);
    for (int j = 0; j <= N; ++j)
      for (int i = 0; i <= N; ++i) {
        Vector3s pos = origin + Vector3s(dx * i, dx * j, 0.0);
        AddVertex(i, j, k, pos, Vector3s(-b(0), -b(1), b(2) * sign[kidx]));

        pos = origin + Vector3s(dx * i + 2.0 * b(0) + radius, dx * j, 0.0);
        AddVertex(i + N + 1, j, k, pos,
                  Vector3s(b(0), -b(1), b(2) * sign[kidx]));

        pos = origin + Vector3s(dx * i + 2.0 * b(0) + radius,
                                dx * j + 2.0 * b(1) + radius, 0.0);
        AddVertex(i + N + 1, j + N + 1, k, pos,
                  Vector3s(b(0), b(1), b(2) * sign[kidx]));

        pos = origin + Vector3s(dx * i, dx * j + 2.0 * b(1) + radius, 0.0);
        AddVertex(i, j + N + 1, k, pos,
                  Vector3s(-b(0), b(1), b(2) * sign[kidx]));
      }
    // corners
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        AddFace(translateIndices(i, j, k), translateIndices(i + 1, j + 1, k),
                translateIndices(i, j + 1, k), kidx == 0);
        AddFace(translateIndices(i, j, k), translateIndices(i + 1, j, k),
                translateIndices(i + 1, j + 1, k), kidx == 0);

        AddFace(translateIndices(i, j + N + 1, k),
                translateIndices(i + 1, j + N + 2, k),
                translateIndices(i, j + N + 2, k), kidx == 0);
        AddFace(translateIndices(i, j + N + 1, k),
                translateIndices(i + 1, j + N + 1, k),
                translateIndices(i + 1, j + N + 2, k), kidx == 0);

        AddFace(translateIndices(i + N + 1, j + N + 1, k),
                translateIndices(i + N + 2, j + N + 2, k),
                translateIndices(i + N + 1, j + N + 2, k), kidx == 0);
        AddFace(translateIndices(i + N + 1, j + N + 1, k),
                translateIndices(i + N + 2, j + N + 1, k),
                translateIndices(i + N + 2, j + N + 2, k), kidx == 0);

        AddFace(translateIndices(i + N + 1, j, k),
                translateIndices(i + N + 2, j + 1, k),
                translateIndices(i + N + 1, j + 1, k), kidx == 0);
        AddFace(translateIndices(i + N + 1, j, k),
                translateIndices(i + N + 2, j, k),
                translateIndices(i + N + 2, j + 1, k), kidx == 0);
      }
    // sides
    for (int i = 0; i < N; ++i) {
      AddFace(translateIndices(i, N, k), translateIndices(i + 1, N + 1, k),
              translateIndices(i, N + 1, k), kidx == 0);
      AddFace(translateIndices(i, N, k), translateIndices(i + 1, N, k),
              translateIndices(i + 1, N + 1, k), kidx == 0);

      AddFace(translateIndices(N, i, k), translateIndices(N + 1, i + 1, k),
              translateIndices(N, i + 1, k), kidx == 0);
      AddFace(translateIndices(N, i, k), translateIndices(N + 1, i, k),
              translateIndices(N + 1, i + 1, k), kidx == 0);

      AddFace(translateIndices(i + N + 1, N, k),
              translateIndices(i + N + 2, N + 1, k),
              translateIndices(i + N + 1, N + 1, k), kidx == 0);
      AddFace(translateIndices(i + N + 1, N, k),
              translateIndices(i + N + 2, N, k),
              translateIndices(i + N + 2, N + 1, k), kidx == 0);

      AddFace(translateIndices(N, i + N + 1, k),
              translateIndices(N + 1, i + N + 2, k),
              translateIndices(N, i + N + 2, k), kidx == 0);
      AddFace(translateIndices(N, i + N + 1, k),
              translateIndices(N + 1, i + N + 1, k),
              translateIndices(N + 1, i + N + 2, k), kidx == 0);
    }
    // central
    AddFace(translateIndices(N, N, k), translateIndices(N + 1, N + 1, k),
            translateIndices(N, N + 1, k), kidx == 0);
    AddFace(translateIndices(N, N, k), translateIndices(N + 1, N, k),
            translateIndices(N + 1, N + 1, k), kidx == 0);
  }

  // xz-planes
  for (int kidx = 0; kidx < 2; ++kidx) {
    const int k = ks[kidx];
    Vector3s origin =
        Vector3s(-b(0) - radius, (b(1) + radius) * sign[kidx], -b(2) - radius);
    for (int j = 0; j <= N; ++j)
      for (int i = 0; i <= N; ++i) {
        Vector3s pos = origin + Vector3s(dx * i, 0.0, dx * j);
        AddVertex(i, k, j, pos, Vector3s(-b(0), b(1) * sign[kidx], -b(2)));

        pos = origin + Vector3s(dx * i + 2.0 * b(0) + radius, 0.0, dx * j);
        AddVertex(i + N + 1, k, j, pos,
                  Vector3s(b(0), b(1) * sign[kidx], -b(2)));

        pos = origin + Vector3s(dx * i + 2.0 * b(0) + radius, 0.0,
                                dx * j + 2.0 * b(2) + radius);
        AddVertex(i + N + 1, k, j + N + 1, pos,
                  Vector3s(b(0), b(1) * sign[kidx], b(2)));

        pos = origin + Vector3s(dx * i, 0.0, dx * j + 2.0 * b(2) + radius);
        AddVertex(i, k, j + N + 1, pos,
                  Vector3s(-b(0), b(1) * sign[kidx], b(2)));
      }
    // corners
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        AddFace(translateIndices(i, k, j), translateIndices(i + 1, k, j + 1),
                translateIndices(i, k, j + 1), kidx == 1);
        AddFace(translateIndices(i, k, j), translateIndices(i + 1, k, j),
                translateIndices(i + 1, k, j + 1), kidx == 1);

        AddFace(translateIndices(i, k, j + N + 1),
                translateIndices(i + 1, k, j + N + 2),
                translateIndices(i, k, j + N + 2), kidx == 1);
        AddFace(translateIndices(i, k, j + N + 1),
                translateIndices(i + 1, k, j + N + 1),
                translateIndices(i + 1, k, j + N + 2), kidx == 1);

        AddFace(translateIndices(i + N + 1, k, j + N + 1),
                translateIndices(i + N + 2, k, j + N + 2),
                translateIndices(i + N + 1, k, j + N + 2), kidx == 1);
        AddFace(translateIndices(i + N + 1, k, j + N + 1),
                translateIndices(i + N + 2, k, j + N + 1),
                translateIndices(i + N + 2, k, j + N + 2), kidx == 1);

        AddFace(translateIndices(i + N + 1, k, j),
                translateIndices(i + N + 2, k, j + 1),
                translateIndices(i + N + 1, k, j + 1), kidx == 1);
        AddFace(translateIndices(i + N + 1, k, j),
                translateIndices(i + N + 2, k, j),
                translateIndices(i + N + 2, k, j + 1), kidx == 1);
      }
    // sides
    for (int i = 0; i < N; ++i) {
      AddFace(translateIndices(i, k, N), translateIndices(i + 1, k, N + 1),
              translateIndices(i, k, N + 1), kidx == 1);
      AddFace(translateIndices(i, k, N), translateIndices(i + 1, k, N),
              translateIndices(i + 1, k, N + 1), kidx == 1);

      AddFace(translateIndices(N, k, i), translateIndices(N + 1, k, i + 1),
              translateIndices(N, k, i + 1), kidx == 1);
      AddFace(translateIndices(N, k, i), translateIndices(N + 1, k, i),
              translateIndices(N + 1, k, i + 1), kidx == 1);

      AddFace(translateIndices(i + N + 1, k, N),
              translateIndices(i + N + 2, k, N + 1),
              translateIndices(i + N + 1, k, N + 1), kidx == 1);
      AddFace(translateIndices(i + N + 1, k, N),
              translateIndices(i + N + 2, k, N),
              translateIndices(i + N + 2, k, N + 1), kidx == 1);

      AddFace(translateIndices(N, k, i + N + 1),
              translateIndices(N + 1, k, i + N + 2),
              translateIndices(N, k, i + N + 2), kidx == 1);
      AddFace(translateIndices(N, k, i + N + 1),
              translateIndices(N + 1, k, i + N + 1),
              translateIndices(N + 1, k, i + N + 2), kidx == 1);
    }
    // central
    AddFace(translateIndices(N, k, N), translateIndices(N + 1, k, N + 1),
            translateIndices(N, k, N + 1), kidx == 1);
    AddFace(translateIndices(N, k, N), translateIndices(N + 1, k, N),
            translateIndices(N + 1, k, N + 1), kidx == 1);
  }

  // yz-planes
  for (int kidx = 0; kidx < 2; ++kidx) {
    const int k = ks[kidx];
    Vector3s origin =
        Vector3s((b(0) + radius) * sign[kidx], -b(1) - radius, -b(2) - radius);
    for (int j = 0; j <= N; ++j)
      for (int i = 0; i <= N; ++i) {
        Vector3s pos = origin + Vector3s(0.0, dx * i, dx * j);
        AddVertex(k, i, j, pos, Vector3s(b(0) * sign[kidx], -b(1), -b(2)));

        pos = origin + Vector3s(0.0, dx * i + 2.0 * b(1) + radius, dx * j);
        AddVertex(k, i + N + 1, j, pos,
                  Vector3s(b(0) * sign[kidx], b(1), -b(2)));

        pos = origin + Vector3s(0.0, dx * i + 2.0 * b(1) + radius,
                                dx * j + 2.0 * b(2) + radius);
        AddVertex(k, i + N + 1, j + N + 1, pos,
                  Vector3s(b(0) * sign[kidx], b(1), b(2)));

        pos = origin + Vector3s(0.0, dx * i, dx * j + 2.0 * b(2) + radius);
        AddVertex(k, i, j + N + 1, pos,
                  Vector3s(b(0) * sign[kidx], -b(1), b(2)));
      }
    // corners
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        AddFace(translateIndices(k, i, j), translateIndices(k, i + 1, j + 1),
                translateIndices(k, i, j + 1), kidx == 0);
        AddFace(translateIndices(k, i, j), translateIndices(k, i + 1, j),
                translateIndices(k, i + 1, j + 1), kidx == 0);

        AddFace(translateIndices(k, i, j + N + 1),
                translateIndices(k, i + 1, j + N + 2),
                translateIndices(k, i, j + N + 2), kidx == 0);
        AddFace(translateIndices(k, i, j + N + 1),
                translateIndices(k, i + 1, j + N + 1),
                translateIndices(k, i + 1, j + N + 2), kidx == 0);

        AddFace(translateIndices(k, i + N + 1, j + N + 1),
                translateIndices(k, i + N + 2, j + N + 2),
                translateIndices(k, i + N + 1, j + N + 2), kidx == 0);
        AddFace(translateIndices(k, i + N + 1, j + N + 1),
                translateIndices(k, i + N + 2, j + N + 1),
                translateIndices(k, i + N + 2, j + N + 2), kidx == 0);

        AddFace(translateIndices(k, i + N + 1, j),
                translateIndices(k, i + N + 2, j + 1),
                translateIndices(k, i + N + 1, j + 1), kidx == 0);
        AddFace(translateIndices(k, i + N + 1, j),
                translateIndices(k, i + N + 2, j),
                translateIndices(k, i + N + 2, j + 1), kidx == 0);
      }
    // sides
    for (int i = 0; i < N; ++i) {
      AddFace(translateIndices(k, i, N), translateIndices(k, i + 1, N + 1),
              translateIndices(k, i, N + 1), kidx == 0);
      AddFace(translateIndices(k, i, N), translateIndices(k, i + 1, N),
              translateIndices(k, i + 1, N + 1), kidx == 0);

      AddFace(translateIndices(k, N, i), translateIndices(k, N + 1, i + 1),
              translateIndices(k, N, i + 1), kidx == 0);
      AddFace(translateIndices(k, N, i), translateIndices(k, N + 1, i),
              translateIndices(k, N + 1, i + 1), kidx == 0);

      AddFace(translateIndices(k, i + N + 1, N),
              translateIndices(k, i + N + 2, N + 1),
              translateIndices(k, i + N + 1, N + 1), kidx == 0);
      AddFace(translateIndices(k, i + N + 1, N),
              translateIndices(k, i + N + 2, N),
              translateIndices(k, i + N + 2, N + 1), kidx == 0);

      AddFace(translateIndices(k, N, i + N + 1),
              translateIndices(k, N + 1, i + N + 2),
              translateIndices(k, N, i + N + 2), kidx == 0);
      AddFace(translateIndices(k, N, i + N + 1),
              translateIndices(k, N + 1, i + N + 1),
              translateIndices(k, N + 1, i + N + 2), kidx == 0);
    }
    // central
    AddFace(translateIndices(k, N, N), translateIndices(k, N + 1, N + 1),
            translateIndices(k, N, N + 1), kidx == 0);
    AddFace(translateIndices(k, N, N), translateIndices(k, N + 1, N),
            translateIndices(k, N + 1, N + 1), kidx == 0);
  }
}
