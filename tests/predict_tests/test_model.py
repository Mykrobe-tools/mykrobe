import sys
import json

from unittest import TestCase

sys.path.append(".")

from mykrobe.predict import MykrobePredictorSusceptibilityResult


class MykrobePredictResultsTest(TestCase):

    def setUp(self):
        pass

    def teardown(self):
        pass

    def test_document_from_json(self):
    	temp_json = {'susceptibility':{'Rifampicin':{"predict" : 'R'}}}
    	result = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json))
    	assert result.susceptibility["Rifampicin"] == {"predict" : 'R'}

    	result2 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json))
    	assert result is not result2
    	assert result == result2


class MykrobePredictResultsTest(TestCase):

    def setUp(self):
    	temp_json = {'susceptibility':{'Rifampicin':{"predict" : 'R'}}}
    	self.result = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json))
    	temp_json2 = {'susceptibility':{'Rifampicin':{"predict" : 'R'}}}
    	self.result2 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json2))
    	temp_json3 = {'susceptibility':{'Rifampicin':{"predict" : 'S'}}}
    	self.result3 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json3))
    	temp_json4 = {'susceptibility':{'Quin':{"predict" : 'S'}}}
    	self.result4 = MykrobePredictorSusceptibilityResult.from_json(json.dumps(temp_json4))

    def teardown(self):
        pass

    def test_document_diff(self):
    	assert self.result.diff(self.result2) == {}
    	assert self.result.diff(self.result3) == {"Rifampicin" : {"predict" :("R", "S")}}

    def test_document_diff2(self):
    	assert self.result.diff(self.result4) == {'Rifampicin': {'predict': ('R', 'NA')}, 'Quin': {'predict': ('NA', 'S')}}

