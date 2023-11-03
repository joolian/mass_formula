"""
Classes to retrieve the formula of a fragment give an m/z and the possible elements.

If you use it, do not forget to cite the authors from the web API:
ChemCalc: a building block for tomorrow's chemical infrastructure.
Patiny, Luc; Borel, Alain Journal of Chemical Information and Modeling 2013.
DOI:10.1021/ci300563h
"""

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from requests.exceptions import RetryError, ConnectionError, JSONDecodeError


class ChemCalcFormulaFinder:
    """
    Uses the www.chemcalc.org API to find the possible molecular formulas of a mass or masses.
    """
    url = 'https://www.chemcalc.org/chemcalc/em'

    df_columns = [
        'em', 'mf', 'unsat', 'jcampURL', 'error', 'ppm', 'info', 'minMass', 'mfRange', 'numberOfResultsOnly', 'typedResult',
        'minUnsaturation', 'maxUnsaturation', 'useUnsaturation', 'jcampBaseURL', 'monoisotopicMass', 'jcampLink',
        'integerUnsaturation', 'maxMass', 'referenceVersion', 'massRange', 'charge', 'number_results', 'search_number', 'accuracy',
        'error_description'
    ]

    def __init__(
            self,
            typed_result=False,
            use_unsaturation=False,
            min_unsaturation=0,
            max_unsaturation=50,
            integer_unsaturation=False,
            jcamp_link=True
    ):
        self._typed_result = typed_result
        self._use_unsaturation = use_unsaturation
        self._min_unsaturation = min_unsaturation
        self._max_unsaturation = max_unsaturation
        self._integer_unsaturation = integer_unsaturation
        self._jcamp_link = jcamp_link
        self._session = requests.Session()
        adapter = HTTPAdapter(
            max_retries=Retry(
                total=4, backoff_factor=1,
                allowed_methods=None, status_forcelist=[429, 500, 502, 503, 504])
        )
        self._session.mount("http://", adapter)

    def results_meta_data(self, json):
        meta = {key: value for (key, value) in json.items() if not isinstance(value, (dict, list))}
        return {**json['options'], **meta, **json['error']}

    def json_to_dataframe(self, formulas_json):
        """ converts the json returned by the formula search to a pandas DataFrame"""
        df = pd.DataFrame(columns=ChemCalcFormulaFinder.df_columns)
        if formulas_json['error']['error_description'] is not None:
            del formulas_json['error']['fatal']
            return pd.concat([
                df,
                pd.DataFrame([formulas_json['error']])
            ])

        meta_data = self.results_meta_data(formulas_json)
        if not formulas_json['results']:
            return pd.concat([df, pd.DataFrame(meta_data)])
        else:
            df = pd.concat([
                df,
                pd.DataFrame.from_dict(formulas_json['results'])
            ])
            df[list(meta_data.keys())] = list(meta_data.values())
            return df

    def get_formulas(self, mass, accuracy, mf_range, charge):
        """finds the possible formulas for a mass"""
        mf_range = f'{mf_range}({charge})'
        params = {
            'mfRange': mf_range,
            'numberOfResultsOnly': False,
            'typedResult': self._typed_result,
            'useUnsaturation': self._use_unsaturation,
            'minUnsaturation': self._min_unsaturation,
            'maxUnsaturation': self._max_unsaturation,
            'jcampBaseURL': 'http://www.chemcalc.org/service/jcamp/',
            'monoisotopicMass': mass,
            'jcampLink': self._jcamp_link,
            'integerUnsaturation': self._integer_unsaturation,
            'referenceVersion': '2013',
            'massRange': accuracy
        }
        try:
            results = self._session.get(ChemCalcFormulaFinder.url, params=params, timeout=5).json()
            results['error'] = {'error_description': None}
            return results
        except (RetryError, ConnectionError) as err:
            fatal_error = True
            error_description = repr(err)
        except JSONDecodeError as err:
            fatal_error = False
            error_description = repr(err)
        return {
            'error': {
                'fatal': fatal_error,
                'error_description': error_description,
                'em': mass,
                'accuracy': accuracy,
                'mfRange': mf_range,
                'charge': charge
            }
        }

    def formulas(self, mass, accuracy, mf_range, charge):
        """

        :param mass: scalar or array of masses
        :param accuracy: scalar or array of mass accuracies expressed as absolute mass
        :param mf_range: String or array of strings listing the number and symbol of each isotope to be used.
        :param charge: Scalar or array of charges on the mass or masses
        :return: pandas DataFrame
        """
        try:
            params = zip(mass, accuracy, mf_range, charge, strict=True)
        except TypeError:
            params = zip([mass], [accuracy], [mf_range], [charge])
        results = []
        for index, (mz, accuracy, mf_range, charge) in enumerate(params):
            formulas_json = self.get_formulas(mz, accuracy, mf_range, charge)
            result = self.json_to_dataframe(formulas_json)
            result['search_number'] = index
            results.append(result)
            if formulas_json['error'].get('fatal_error', False):
                break
        return pd.concat(results)


def ppm_error(mass, mass_error):
    """Calculate the ppm mass error of a mass given the mass and the absolute mass error"""
    return mass_error / mass * 1e6


def ppm_to_mass(mass, ppm_error):
    """Calculates the absolute mass error of a mass given the mass and the ppm mass error"""
    return mass * ppm_error / 1e6


if __name__ == '__main__':
    trichoromethane_monoisotopic_mass = 117.914383
    electron_mass = 9.1093837e-28
    avogadro_constant = 6.02214076e23
    electron_mass = avogadro_constant * electron_mass
    print(f'Electron mass: {electron_mass}')
    print(f'Cl3CH: {trichoromethane_monoisotopic_mass - electron_mass}')
    mz = trichoromethane_monoisotopic_mass - electron_mass

    finder = ChemCalcFormulaFinder(use_unsaturation=True, integer_unsaturation=True)
    mass_accuracy_ppm = 5
    mass_accuracy = ppm_to_mass(mz, mass_accuracy_ppm)
    mf_range = 'C0-100H0-202N0-10O0-10F0-3Cl0-3Br0-1'
    charge = '+'
    result = finder.formulas(mz, mass_accuracy, mf_range, charge)
    print(result)

    masses = [50, 150]
    mass_accuracy = [mass_accuracy, mass_accuracy]
    mf_range = [mf_range, '']
    charge = ['+', '+']

    data = pd.DataFrame(data={'mz': masses, 'mass_accuracy': mass_accuracy, 'mf_range': mf_range, 'charge': charge})
    result = finder.formulas(data.mz, data.mass_accuracy, data.mf_range, data.charge)
    result.to_csv('results.csv')
    print(result)
