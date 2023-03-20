# import click
# from click.testing import CliRunner
#
# @click.command()
# @click.option('-w', '--width', type=int, default=0)
# @click.option('--option2')
# @click.argument('argument')
# def app(width, option2, argument):
#     click.echo("params: {} {} {}".format(width, option2, argument))
#
#
# def test_app():
#     # Shows an example of how to setup CliRunner(), invoke a command line app, and process the output
#     runner = CliRunner()
#     result = runner.invoke(app, ["-w", 3, "--option2", 4, "arg"])
#     assert result.exit_code == 0
#     assert result.output == 'params: 3 4 arg\n'
#     print ('done')
#
