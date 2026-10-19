# GitHub代码包使用说明

把本文件所在目录中的内容上传至GitHub，目录根部应直接看到 `README.md`、`run_analysis.py`、`requirements.txt`、`analysis/`、`audit/` 和 `scripts/`。保留隐藏的 `.gitignore`。

本包包括独立原始数据重建及完整病例模型、原始seed 20261011多重插补代码、seed 20261017/20261018独立核验、全精度汇总结果、森林图代码、输入文件清单和依赖说明。不是仅凭汇总表画图的代码包。

不附带原始个体数据；官方来源、四份准确文件名、ZIP内部文件名和SHA-256清单详见英文README。数据下载的 `.tab` 和 `.dta` 格式不可仅靠改扩展名互换。`clhls_1998_2018_longitudinal_dataset_released_version1.dta` 不是本代码实际使用的2018横断面输入。

最简运行：

```bash
python -m pip install -r requirements.txt
python run_analysis.py verify-archive
python run_analysis.py figure
```

以上两项不用原始数据，核验归档结果并生成图。若要重新估计，把四份授权原始文件放入 `inputs/` 后运行：

```bash
python run_analysis.py main
```

原始多重插补需额外安装statsmodels：

```bash
python -m pip install -r requirements-original-mi.txt
python run_analysis.py mi-original
python run_analysis.py mi-audit
```

`all` 可依次执行全部步骤。原始seed 20261011的statsmodels确切版本没有保存在找回的审计元数据中，因此没有虚构完整历史锁定环境；独立重建所用NumPy/pandas/SciPy版本已固定。

这次打包实际通过了语法检查、38行归档效应值/P值核验、四个主分析图点核验、数值求解与结局缺失逻辑检查，并成功生成森林图。当前工作区没有四份原始数据，因此本轮没有新跑完整原始数据或MI模型；此前独立复核的结果与代码均保留。最终公开可复现版本建议由持有原始数据的作者完整运行并留存日志。

图1使用五项ADL主分析：HAALSI 1199/105、CLHLS 3824/623。排除两例原始回忆回答不知道后的1197/104只用于敏感性分析。不要按P值选择版本。

原始数据及运行生成的任何 `.pkl` 个体数据检查点均不要上传。仓库尚未代您公开，也没有代作者选择开源许可证。取得真实GitHub链接及Zenodo DOI后再填入论文。
