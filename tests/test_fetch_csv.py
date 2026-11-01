import pytest
from occulta import survey


class R:
    def __init__(self, text): self.text = text


def setup(monkeypatch, tmp_path, text):
    monkeypatch.setattr(survey, 'CACHE', tmp_path)
    monkeypatch.setattr(survey, 'LOGS', tmp_path)
    monkeypatch.setattr(survey.time, 'sleep', lambda s: None)
    monkeypatch.setattr(survey.requests, 'get', lambda *a, **k: R(text))


def test_missing_column_uses_cache(monkeypatch, tmp_path):
    setup(monkeypatch, tmp_path, '<html>rate limited</html>')
    (tmp_path / 'toi.csv').write_text('TIC ID,TOI\n1,2\n')
    assert 'TIC ID' in survey.fetch_csv('toi', 'u', 'TIC ID')
    assert (tmp_path / 'toi.csv').read_text() == 'TIC ID,TOI\n1,2\n'


def test_missing_column_no_cache_raises(monkeypatch, tmp_path):
    setup(monkeypatch, tmp_path, '<html>rate limited</html>')
    with pytest.raises(RuntimeError, match='TIC ID'):
        survey.fetch_csv('toi', 'u', 'TIC ID')
    assert not (tmp_path / 'toi.csv').exists()
