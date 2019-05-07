import argparse
import numpy as np
from argparse import RawTextHelpFormatter


def get_heron_area(a, b, c):

    x = np.linalg.norm((b - a), 2)
    y = np.linalg.norm((c - a), 2)
    z = np.linalg.norm((c - b), 2)
    s = (x + y + z) * 0.5

    return (s * (s - x) * (s - y) * (s - z)) ** 0.5


def calc_A_mixed(vertices, triangles):

    numv = vertices.shape[0]
    numt = triangles.shape[0]

    A_mixed = np.zeros((numv, numt))

    mean_curvature_normal_operator = np.zeros((numv, numt, 3))

    for i in range(numv):

        req_t = triangles[(triangles[:, 0] == i) | (
            triangles[:, 1] == i) | (triangles[:, 2] == i)]

        for j in range(len(req_t)):

            tid = np.where(np.all(triangles == req_t[j], axis=1))

            nbhr = [v for v in req_t[j] if v != i]

            vec1 = (vertices[nbhr[0]] - vertices[i]) / \
                np.linalg.norm(vertices[nbhr[0]] - vertices[i], 2)
            vec2 = (vertices[nbhr[1]] - vertices[i]) / \
                np.linalg.norm(vertices[nbhr[1]] - vertices[i], 2)
            angle_at_x = np.arccos(np.dot(vec1, vec2))

            if angle_at_x > np.pi / 2:
                A_mixed[i, tid] = get_heron_area(
                    vertices[i], vertices[nbhr[0]], vertices[nbhr[1]]) / 2
                continue

            vec1a = (vertices[i] - vertices[nbhr[0]]) / \
                np.linalg.norm(vertices[i] - vertices[nbhr[0]], 2)
            vec2a = (vertices[nbhr[1]] - vertices[nbhr[0]]) / \
                np.linalg.norm(vertices[nbhr[1]] - vertices[nbhr[0]], 2)

            inner_prod = np.dot(vec1a, vec2a)
            angle1 = np.arccos(inner_prod)

            if angle1 > np.pi / 2:
                A_mixed[i, tid] = get_heron_area(
                    vertices[i], vertices[nbhr[0]], vertices[nbhr[1]]) / 4
                continue

            vec1b = (vertices[i] - vertices[nbhr[1]]) / \
                np.linalg.norm(vertices[i] - vertices[nbhr[1]], 2)
            vec2b = (vertices[nbhr[0]] - vertices[nbhr[1]]) / \
                np.linalg.norm(vertices[nbhr[0]] - vertices[nbhr[1]], 2)

            inner_prod = np.dot(vec1b, vec2b)
            angle2 = np.arccos(inner_prod)

            if angle2 > np.pi / 2:
                A_mixed[i, tid] = get_heron_area(
                    vertices[i], vertices[nbhr[0]], vertices[nbhr[1]]) / 4
                continue

            cot_1 = 1 / np.tan(angle1)
            cot_2 = 1 / np.tan(angle2)

            A_v_of_tid = 0.125 * ((cot_1 * np.linalg.norm(vertices[i] - vertices[nbhr[
                1]], 2)**2) + (cot_2 * np.linalg.norm(vertices[i] - vertices[nbhr[0]], 2)**2))

            mean_curvature_normal_operator_at_v_t = ((1 / np.tan(angle1)) * (
                vertices[i] - vertices[nbhr[1]])) + ((1 / np.tan(angle2)) * (vertices[i] - vertices[nbhr[0]]))

            A_mixed[i, tid] = A_v_of_tid
            mean_curvature_normal_operator[
                i, tid] = mean_curvature_normal_operator_at_v_t

    A_mixed = np.sum(A_mixed, axis=1)
    # Set zeros in A_mixed to very small values
    A_mixed[A_mixed == 0] = 10 ** -40
    mean_curvature_normal_operator = (
        (1 / (2 * A_mixed)) * np.sum(mean_curvature_normal_operator, axis=1).T).T

    return A_mixed, mean_curvature_normal_operator


def get_mean_curvature(mean_curvature_normal_operator_vector):
    K_H = 0.5 * \
        np.linalg.norm(mean_curvature_normal_operator_vector, 2, axis=1)
    return K_H


def get_gaussian_curvature(vertices, triangles, A_mixed):
    numv = vertices.shape[0]
    numt = triangles.shape[0]
    K_G = np.zeros(numv)
    for i in range(numv):
        sum_theta = 0
        req_t = triangles[(triangles[:, 0] == i) | (
            triangles[:, 1] == i) | (triangles[:, 2] == i)]

        for j in range(req_t.shape[0]):

            nbhrs = [v for v in req_t[j] if v != i]
            vec1 = vertices[nbhrs[0]] - vertices[i]
            vec1 = vec1 / np.linalg.norm(vec1, 2)
            vec2 = vertices[nbhrs[1]] - vertices[i]
            vec2 = vec2 / np.linalg.norm(vec2, 2)
            angle = np.arccos(np.dot(vec1, vec2))
            sum_theta += angle

        K_G[i] = ((2 * np.pi) - sum_theta) / A_mixed[i]
    return K_G


def get_principal_curvatures(K_H, K_G):
    numv = vertices.shape[0]
    numt = triangles.shape[0]
    zeros = np.zeros(numv)
    delx = np.sqrt(np.max(np.vstack((K_H**2 - K_G, zeros)), axis=0))
    K_1 = K_H + delx
    K_2 = K_H - delx
    return K_1, K_2


def read_off(file):
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')

    n_verts, n_faces, n_dontknow = tuple(
        [int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')]
             for i_vert in range(n_verts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:]
             for i_face in range(n_faces)]
    return np.array(verts, dtype=np.float64), np.array(faces, dtype=np.int64)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="First, use '--mode 0' to generate files for containing value of the \noperator and then plot the operatore using '--mode 1' and '--mesh'. For \n--ops, the operations are encoded as: \n\t- 1: Mean Curvature\n\t- 2: Gaussian Curvature\n\t- 3: Principal Curvatures\nNote that each operation is performed for all the input files.", formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "--mode", help="specifies mode for program:\n- 0: For computation\n- 1: For Plotting", type=int, required=True, default=0)
    parser.add_argument("i", help="path to input file(s)", type=str, nargs="+")
    parser.add_argument(
        "--ops", help="number to denote all the operations to be \nperformed on each file.", type=int)
    parser.add_argument(
        "--mesh", help="mesh on which curvatures were calculated (redundant in computation mode)", default=None)
    parser.add_argument(
        "--save", help="flag for saving the plot (redundant in computation mode)", default=False, action="store_true")
    parser.add_argument(
        "--title", help="title of plot (redundant in computation mode)", default="My Plot")

    args = parser.parse_args()
    if args.mode == 0:

        for inp in args.i:

            ops = [int(d) for d in str(args.ops)]
            if ops is None:
                parser.error("--mode 0 requires --ops")
            mesh_file = inp
            f = open(mesh_file)
            vertices, triangles = read_off(f)

            A_mixed = None
            mean_curvature_normal_operator_vec = None

            A_mixed, mean_curvature_normal_operator_vec = calc_A_mixed(
                vertices, triangles)

            K_H = None
            K_G = None
            K_1 = None
            K_2 = None

            for op in ops:
                if op == 1 and K_H is None:
                    K_H = get_mean_curvature(
                        mean_curvature_normal_operator_vec)
                    np.save("./" + inp.split(".")[0] + "_KH.npy", K_H)
                    print("[DiffGeoOps]: Mean Curvature for", inp,
                          "saved to", "./" + inp.split(".")[0] + "_KH.npy")

                elif op == 2 and K_G is None:
                    K_G = get_gaussian_curvature(vertices, triangles, A_mixed)
                    np.save("./" + inp.split(".")[0] + "_KG.npy", K_G)
                    print("[DiffGeoOps]: Gaussian Curvature for", inp,
                          "saved to", "./" + inp.split(".")[0] + "_KG.npy")

                elif op == 3:
                    if K_H is None:
                        K_H = get_mean_curvature(
                            mean_curvature_normal_operator_vec)
                    if K_G is None:
                        K_G = get_gaussian_curvature(
                            vertices, triangles, A_mixed)
                    K_1, K_2 = get_principal_curvatures(K_H, K_G)
                    np.save("./" + inp.split(".")[0] + "_K1.npy", K_1)
                    np.save("./" + inp.split(".")[0] + "_K2.npy", K_2)
                    print("[DiffGeoOps]: Principal Curvature 1 for", inp,
                          "saved to", "./" + inp.split(".")[0] + "_K1.npy")
                    print("[DiffGeoOps]: Principal Curvature 2 for", inp,
                          "saved to", "./" + inp.split(".")[0] + "_K2.npy")

    elif args.mode == 1:

        if len(args.i) > 1:
            parser.error("Multiple inputs only allowed in computation mode!")

        from mayavi import mlab

        if args.i[0].split(".")[1] != "npy":
            raise RunTimeError("Plotting requires .npy files!")

        f = np.load(args.i[0])
        if args.mesh is None:
            parser.error("Plotting requires --mesh")
        mlab.figure(args.title, size=(600, 600))
        mesh_file = open(args.mesh)
        vertices, triangles = read_off(mesh_file)

        x, y, z = vertices[:, 0], vertices[:, 1], vertices[:, 2]

        mesh = mlab.triangular_mesh(x, y, z, triangles,
                                    representation='wireframe', opacity=0)

        mesh.mlab_source.dataset.point_data.scalars = f
        mesh.mlab_source.dataset.point_data.scalars.name = 'Point data'

        mesh.mlab_source.update()
        mesh.parent.update()

        mesh2 = mlab.pipeline.set_active_attribute(mesh,
                                                   point_scalars='Point data')
        s2 = mlab.pipeline.surface(mesh2)
        s2.actor.mapper.interpolate_scalars_before_mapping = True
        mlab.colorbar(s2, title='Curvature\n', orientation='vertical')

        if args.save:
            mlab.savefig(args.i[0].split(".")[0] + ".png")
            print("[DiffGeoOps]: Plot for", args.i[0],
                  "saved to", args.i[0].split(".")[0] + ".png")
        print("[DiffGeoOps]: Showing plot for", args.i[0])
        mlab.show()

    else:
        parser.error("Flag not recognized. Please use -h for usage.")
