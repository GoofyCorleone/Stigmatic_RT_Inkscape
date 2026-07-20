from pytest import approx

from inkscape_raytracing.raytracing.geometry.cubic_bezier import cubic_real_roots


def test_cubic_real_roots():
    def roots_sorted(a0, a1, a2, a3):
        return sorted(cubic_real_roots(a0, a1, a2, a3))

    # true cubic
    assert roots_sorted(-12, 22, -12, 2) == approx([1, 2, 3])
    assert roots_sorted(-0, 1, -2, 1) == approx([0, 1])
    assert roots_sorted(-8, 0, 0, 1) == approx([2])
    assert roots_sorted(1, 2, 0, 1) == approx([-0.453398], abs=1e-5)
    assert roots_sorted(0, 0, 0, 1) == approx([0])
    assert roots_sorted(-1, 3, -3, 1) == approx([1])
    # quadratic
    assert roots_sorted(1, -2, 1, 0) == approx([1])
    assert roots_sorted(-1, 0, 1, 0) == approx([-1, 1])
    assert roots_sorted(1, 0, 1, 0) == []
    # linear
    assert roots_sorted(1, 2, 0, 0) == approx([-0.5])
    assert roots_sorted(1, 0, 0, 0) == []
    assert roots_sorted(0, 0, 0, 0) == []
