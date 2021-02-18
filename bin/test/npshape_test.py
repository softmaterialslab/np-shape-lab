import testbook
import os
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

