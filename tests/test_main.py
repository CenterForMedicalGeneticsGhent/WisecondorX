import unittest
from typer.testing import CliRunner
from wisecondorx.main import app, VERSION

runner = CliRunner()


class TestMain(unittest.TestCase):
    def test_version_flag(self):
        result = runner.invoke(app, ["--version"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn(f"WisecondorX version: {VERSION}", result.stdout)

    def test_version_flag_short(self):
        result = runner.invoke(app, ["-v"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn(f"WisecondorX version: {VERSION}", result.stdout)


if __name__ == "__main__":
    unittest.main()
