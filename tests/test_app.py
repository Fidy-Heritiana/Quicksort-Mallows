"""Exercise all pages and check that displayed parameters follow submitted runs."""
import unittest

from streamlit.testing.v1 import AppTest

from explore import ROOT


class ApplicationChecks(unittest.TestCase):
    def test_pages_and_submitted_experiment(self):
        app = AppTest.from_file(str(ROOT / "streamlit_app.py"), default_timeout=30).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="saved_phi").value, "0.5")
        self.assertEqual(app.metric[0].value, "250,658.337")
        app.selectbox(key="saved_phi").set_value("0.8").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.metric[0].value, "102,391.100")

        app.radio(key="page").set_value("Moment convergence").run()
        self.assertFalse(app.exception)
        app.multiselect(key="moment_phis").set_value([]).run()
        self.assertFalse(app.exception)
        self.assertTrue(app.info)
        app.multiselect(key="moment_phis").set_value(["0.5"]).run()
        app.selectbox(key="moment_metric").set_value("sixth_standardised_moment").run()
        self.assertFalse(app.exception)

        app.radio(key="page").set_value("Model and verification").run()
        self.assertFalse(app.exception)

        app.radio(key="page").set_value("New simulation").run()
        app.number_input(key="custom_n").set_value(3)
        app.slider(key="custom_reps").set_value(500)
        app.button[0].click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.metric[0].value, "2.714")
        self.assertEqual(app.session_state["custom_result"]["n"], 3)
        # Unsubmitted controls must not relabel the last completed calculation.
        app.number_input(key="custom_n").set_value(30).run()
        self.assertEqual(app.session_state["custom_result"]["n"], 3)
        self.assertEqual(app.metric[0].value, "2.714")


if __name__ == "__main__":
    unittest.main()
