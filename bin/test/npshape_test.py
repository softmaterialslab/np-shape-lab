import testbook
import os, shutil
import svgwrite
import png

@testbook.testbook('npshape-frontend.ipynb', execute=True)
def test_combine_4_PNGs(tb):
    f = open('test/resources/test.png', 'wb')  # binary mode is important
    w = png.Writer(256, 1, greyscale=True)
    w.write(f, [range(256)])
    f.close()
    if not os.path.exists('snapshot_png'):
        os.mkdir("snapshot_png")

    func = tb.ref("combine_4_PNGs")
    func("test/resources/test.png", "test/resources/test.png", "test/resources/test.png", "test/resources/test.png", "test/resources/output.png")
    assert os.path.exists("test/resources/output.png") == 1
    os.remove("test/resources/output.png")

@testbook.testbook('npshape-frontend.ipynb', execute=True)
def test_convert_SVG_PNG(tb):
    func = tb.ref("convert_SVG_PNG")

    dwg = svgwrite.Drawing('test/resources/test.svg', profile='tiny', size=(270, 270))
    dwg.add(dwg.line((0, 0), (10, 0), stroke=svgwrite.rgb(10, 10, 16, '%')))
    dwg.add(dwg.text('Test', insert=(0, 0.2), fill='red'))
    dwg.save()
    if not os.path.exists('snapshot_png'):
        os.mkdir("snapshot_png")

    func("test/resources/test.svg", "test/resources/output.png")
    assert os.path.exists("test/resources/output.png") == 1
    os.remove("test/resources/output.png")

#@testbook.testbook('npshape-frontend.ipynb', execute=True)
#def test_create_out_dirs(tb):
#    func = tb.ref("createOutDirs")
#    func()
#    assert os.path.exists("outfiles")
#    assert os.path.exists("simulation_snapshots")
#    assert os.path.exists("snapshot_png")
#    shutil.rmtree('outfiles', True)
#    shutil.rmtree('simulation_snapshots', True)
#    shutil.rmtree('snapshot_png', True)

#@testbook.testbook('npshape-frontend.ipynb', execute=True)
#def test_paramerter_append(tb):
#    func = tb.ref("appendParameters")
#    params = 'Maximum Surface Charge(Q in e) is 1\nNumber of patches used: 1\n' \
#             'Patch size (p) is 1\nSurface Charge(Q in e) is 1\nBending Modulus is ' \
#             '1\nStretching Modulus is 1\nSalt Concentration (M) is 1\n'
#    out = func(1, 1, 1, 1, 1, 1, 1, 4)
#    assert params == out




