import json

from mykrobe.predict import MykrobePredictorSusceptibilityResult


def test_document_from_json():
    temp_json = {"susceptibility": {"Rifampicin": {"predict": "R"}}}
    result = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json))
    assert result.susceptibility["Rifampicin"] == {"predict": "R"}

    result2 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json))
    assert result is not result2
    assert result == result2


def test_mykrobe_predict_results():
    temp_json = {"susceptibility": {"Rifampicin": {"predict": "R"}}}
    result = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json))
    temp_json2 = {"susceptibility": {"Rifampicin": {"predict": "R"}}}
    result2 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json2))
    temp_json3 = {"susceptibility": {"Rifampicin": {"predict": "S"}}}
    result3 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json3))
    temp_json4 = {"susceptibility": {"Quin": {"predict": "S"}}}
    result4 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json4))

    assert result.diff(result2) == {}
    assert result.diff(result3) == {"Rifampicin": {"predict": ("R", "S")}}
    assert result.diff(result4) == {
        "Rifampicin": {"predict": ("R", "NA")},
        "Quin": {"predict": ("NA", "S")},
    }
