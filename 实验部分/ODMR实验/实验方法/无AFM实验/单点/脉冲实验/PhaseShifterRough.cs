using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using Controls.Windows;
using HardWares.射频源.Rigol_DSG_3060;
using MathLib.NormalMath.Decimal;
using ODMR_Lab.IO操作;
using ODMR_Lab.ODMR实验;
using ODMR_Lab.基本控件;
using ODMR_Lab.基本控件.一维图表;
using ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM实验;
using ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM实验.单点.脉冲实验;
using ODMR_Lab.实验部分.序列编辑器;
using ODMR_Lab.实验部分.扫描基方法;
using ODMR_Lab.实验部分.扫描基方法.扫描范围;
using ODMR_Lab.设备部分;
using ODMR_Lab.设备部分.射频源_锁相放大器;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using Window = System.Windows.Window;
using Controls.Charts;
using System.Windows.Media;
using ODMR_Lab.实验部分.扫描基方法.扫描任务.多轮一维扫描;
using ODMR_Lab.实验部分.扫描基方法.扫描任务.多轮一维扫描.数据处理方法;
using ODMR_Lab.设备部分.电源;

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM.点实验.脉冲实验
{
    class PhaseShifterRough : PulseExpBase
    {
        public override bool Is1DScanExp { get; set; } = false;
        public override bool Is2DScanExp { get; set; } = false;

        public override string ODMRExperimentName { get; set; } = "相移器角度粗调";

        public override string ODMRExperimentGroupName { get; set; } = "点实验";

        public override string Description { get; set; } = "使用的序列文件名称：PhaseShifter";

        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            //0.Pi脉冲长度(整数),1.T1间隔长度(整数),2.采样循环次数(整数)，3.超时时间（整数）4.微波频率（小数）5.微波功率（小数）
            new Param<int>("电压最小值(V)",0,"VoltageMin"),
            new Param<int>("电压最大值(V)",10,"VoltageMax"),
            new Param<int>("电压点数",20,"VoltagePoints"),
            new Param<int>("测量次数",1000,"LoopCount"),        //外部的循环
            new Param<MultiScanType>("测量循环类型",MultiScanType.正向扫描,"ScanType"),
            new Param<int>("序列循环次数",1000,"SeqLoopCount"),   //板卡内的
            new Param<double>("微波频率(MHz)",2870,"RFFrequency"),
            new Param<double>("微波功率(dBm)",-20,"RFAmplitude"),
            new Param<int>("单点超时时间",10000,"TimeOut"),
            //new Param<bool>("将结果同步到全局脉冲",true,"ToGlobal")
        };
        public override List<ParamB> OutputParams { get; set; } = new List<ParamB>()
        {
        };
        public override List<KeyValuePair<DeviceTypes, Param<string>>> PulseExpDevices { get; set; } = new List<KeyValuePair<DeviceTypes, Param<string>>>()
        {
            new KeyValuePair<DeviceTypes, Param<string>>(DeviceTypes.电源,new Param<string>("相移器电源","","Power"))
        };
        public override List<ChartData1D> D1ChartDatas { get; set; } = new List<ChartData1D>();
        public override List<ChartData2D> D2ChartDatas { get; set; } = new List<ChartData2D>();
        public override List<FittedData1D> D1FitDatas { get; set; } = new List<FittedData1D>();
        protected override List<ODMRExpObject> GetSubExperiments()
        {
            return new List<ODMRExpObject>()
            {
            };
        }

        public override bool IsAFMSubExperiment { get; protected set; } = true;

        public override bool PreConfirmProcedure()
        {
            if (MessageWindow.ShowMessageBox("提示", "是否要继续?此操作将清除原先的实验数据", MessageBoxButton.YesNo, owner: Window.GetWindow(ParentPage)) == MessageBoxResult.Yes)
            {
                return true;
            }
            return false;
        }

        private int CurrentLoop = 0;
        public List<object> ScanEvent(object device, D1NumricScanRangeBase range, double locvalue, int currrentloop, List<Tuple<string, string, double, MultiLoopDataProcessBase>> outputparams, List<object> inputParams)
        {
            //设置locvalue（电压）
            //GlobalPulseParams.SetGlobalPulseLength("RabiTime", (int)locvalue);
            var power = GetDeviceByName("Power") as PowerChannelInfo;
            power.Channel.Voltage = locvalue;

            PulsePhotonPack pack = DoPulseExp("PhaseShifter", GetInputParamValueByName("RFFrequency"), GetInputParamValueByName("RFAmplitude"), GetInputParamValueByName("SeqLoopCount"), 4, GetInputParamValueByName("TimeOut"));

            //光子数
            double signalcount = pack.GetPhotonsAtIndex(0).Sum();
            double refcount = pack.GetPhotonsAtIndex(1).Sum();

            //计算对比度
            double signalcontrast = 1;
            try
            {
                if (refcount != 0)
                    signalcontrast = (signalcount - refcount) / refcount;
            }
            catch (Exception)
            {
            }

            outputparams.Add(new Tuple<string, string, double, MultiLoopDataProcessBase>("相移器电压(V)", "对比度数据", locvalue, new StandardDataProcess()));
            outputparams.Add(new Tuple<string, string, double, MultiLoopDataProcessBase>("信号对比度[(sig-ref)/ref]", "对比度数据", signalcontrast, new StandardDataProcess()));
            outputparams.Add(new Tuple<string, string, double, MultiLoopDataProcessBase>("相移器电压(V)", "荧光数据", locvalue, new StandardDataProcess()));
            outputparams.Add(new Tuple<string, string, double, MultiLoopDataProcessBase>("平均光子数", "荧光数据", refcount, new StandardDataProcess()));
            outputparams.Add(new Tuple<string, string, double, MultiLoopDataProcessBase>("信号光子数", "荧光数据", signalcount, new StandardDataProcess()));
            return new List<object>();
        }

        public override void ODMRExpWithoutAFM()
        {
            MultiScan1DSession<object> Session = new MultiScan1DSession<object>();
            Session.FirstScanEvent = ScanEvent;
            Session.ScanEvent = ScanEvent;
            Session.ScanSource = null;
            Session.PlotEvent = PlotEvent;
            Session.ProgressBarMethod = new Action<object, double>((obj, v) =>
            {
                SetProgress(v);
            });
            Session.SetStateMethod = new Action<object, int, double>((obj, loop, v) =>
            {
                SetExpState("当前扫描轮数:" + loop.ToString() + ",电压值: " + Math.Round(v, 5).ToString());
            });
            D1NumricLinearScanRange range = new D1NumricLinearScanRange(GetInputParamValueByName("VoltageMin"), GetInputParamValueByName("VoltageMax"), GetInputParamValueByName("VoltagePoints"));

            Session.StateJudgeEvent = JudgeThreadEndOrResumeAction;
            Session.BeginScan(GetInputParamValueByName("LoopCount"), GetInputParamValueByName("ScanType"), range, 0, 100);
        }

        private void PlotEvent(List<MultiLoopScanData> list)
        {
            (Get1DChartData("相移器电压(V)", "对比度数据") as NumricChartData1D).Data = MultiLoopScanData.GetAverageData(list, "相移器电压(V)", "对比度数据");
            (Get1DChartData("信号对比度[(sig-ref)/ref]", "对比度数据") as NumricChartData1D).Data = MultiLoopScanData.GetAverageData(list, "信号对比度[(sig-ref)/ref]", "对比度数据");
            
            (Get1DChartData("相移器电压(V)", "荧光数据") as NumricChartData1D).Data = MultiLoopScanData.GetAverageData(list, "相移器电压(V)", "荧光数据");
            (Get1DChartData("平均光子数", "荧光数据") as NumricChartData1D).Data = MultiLoopScanData.GetAverageData(list, "平均光子数", "荧光数据");
            (Get1DChartData("信号光子数", "荧光数据") as NumricChartData1D).Data = MultiLoopScanData.GetAverageData(list, "信号光子数", "荧光数据");

            (Get1DChartData("相移器电压(V)", "方差") as NumricChartData1D).Data = MultiLoopScanData.GetAverageData(list, "相移器电压(V)", "对比度数据");
            (Get1DChartData("平均光子数", "方差") as NumricChartData1D).Data = MultiLoopScanData.GetSigmaData(list, "平均光子数", "荧光数据");
            (Get1DChartData("信号光子数", "方差") as NumricChartData1D).Data = MultiLoopScanData.GetSigmaData(list, "信号光子数", "荧光数据");
            (Get1DChartData("信号对比度[(sig-ref)/ref]", "方差") as NumricChartData1D).Data = MultiLoopScanData.GetSigmaData(list, "信号对比度[(sig-ref)/ref]", "对比度数据");


            UpdatePlotChart();
            UpdatePlotChartFlow(true);
        }

        public override void PreExpEventWithoutAFM()
        {
            //打开微波
            SignalGeneratorChannelInfo RF = GetDeviceByName("RFSource") as SignalGeneratorChannelInfo;
            RF.Device.IsOutOpen = true;

            D1ChartDatas = new List<ChartData1D>()
            {
                new NumricChartData1D("相移器电压(V)","对比度数据",ChartDataType.X),
                new NumricChartData1D("信号对比度[(sig-ref)/ref]","对比度数据",ChartDataType.Y),

                new NumricChartData1D("相移器电压(V)","荧光数据",ChartDataType.X),
                new NumricChartData1D("平均光子数","荧光数据",ChartDataType.Y),
                new NumricChartData1D("信号光子数","荧光数据",ChartDataType.Y),

                new NumricChartData1D("相移器电压(V)","方差",ChartDataType.X),
                new NumricChartData1D("平均光子数","方差",ChartDataType.Y),
                new NumricChartData1D("信号光子数","方差",ChartDataType.Y),
                new NumricChartData1D("信号对比度[(sig-ref)/ref]","方差",ChartDataType.Y),
            };
            UpdatePlotChart();
            Show1DChartData("对比度数据", "相移器电压(V)", "信号对比度[(sig-ref)/ref]");//实验前展示的数据
        }

        //T1拟合函数
        private double RabiFitFunc(double x, double[] ps)
        {
            double a = ps[0];
            double b = ps[1];
            double c = ps[2];
            double d = ps[3];
            return a * Math.Cos(2 * Math.PI / b * (x - c)) + d;
        }

        public override void AfterExpEventWithoutAFM()
        {
            #region 计算通道X 脉冲周期
            var xs = Get1DChartDataSource("相移器电压(V)", "对比度数据");
            var ys_x = Get1DChartDataSource("信号对比度[(sig-ref)/ref]", "对比度数据");

            //设置初值
            double d_x = ys_x.Average();
            double a_x = Math.Abs(ys_x.Min() - ys_x.Max()) / 2;
            double c_x = 10;
            //Pi脉冲时间
            double b_x = 0;
            try
            {
                var fys_x = ys_x.Select(x => x - ys_x.Average()).ToArray();
                Fourier.Forward(fys_x, Enumerable.Repeat(0.0, ys_x.Count).ToArray(), FourierOptions.Matlab);
                fys_x = fys_x.ToList().GetRange(0, (fys_x.Length - 1) / 2).ToArray();
                var freqs_x = Fourier.FrequencyScale(ys_x.Count, (int)1.0 / Math.Abs(xs[0] - xs[1]));
                freqs_x = freqs_x.ToList().GetRange(0, (freqs_x.Length - 1) / 2).ToArray();
                b_x = 1.0 / freqs_x[fys_x.ToList().IndexOf(fys_x.Max())];
                if (double.IsInfinity(b_x)) b_x = 0;
            }
            catch (Exception)
            {
            }

            double[] ps_x = CurveFitting.FitCurveWithFunc(xs, ys_x, new List<double>() { a_x, b_x, c_x, d_x }, new List<double>() { 10, 10, 10, 10, 10 }, RabiFitFunc, AlgorithmType.LevenbergMarquardt, 20000);

            //设置拟合曲线
            var ftxs = new D1NumricLinearScanRange(xs.Min(), xs.Max(), 500).ScanPoints;
            var fitys_x = ftxs.Select(x => RabiFitFunc(x, ps_x)).ToList();
            D1FitDatas.Add(new FittedData1D("a*cos(2*pi/b*(x-c))+d", "x", new List<string>() { "a", "b", "c", "d" }, ps_x.ToList(), "相移器电压(V)", "对比度数据", new NumricDataSeries("拟合曲线X", ftxs, fitys_x) { LineColor = Colors.LightSkyBlue }));
            UpdatePlotChart();
            UpdatePlotChartFlow(true);
            Show1DFittedData("拟合曲线X");

            double b = ps_x[1]; //拟合b，周期
            double c = ps_x[2]; //拟合c，相位
            double vol90 = b/ 4 +c;
            double vol270 = b *3/ 4 + c;
            //处理vol到最小
            vol90 = vol90-b*(int)(vol90 / b);
            vol270 = vol270 - b * (int)(vol270 / b);
            //OutputParams.Add(new Param<double>("0°电压(V)", ps_x[2], "Vol_0"));
            OutputParams.Add(new Param<double>("90°电压(V)", vol90, "Vol_90"));
            OutputParams.Add(new Param<double>("270°电压(V)", vol270, "Vol_270"));

            //计算平均光子计数
            OutputParams.Add(new Param<double>("平均光子计数", Get1DChartDataSource("平均光子数", "荧光数据").Average(), "X_AverageCount"));
            #endregion
            
            
        }

        public override List<ParentPlotDataPack> GetD1PlotPacks()
        {
            List<ParentPlotDataPack> PlotData = new List<ParentPlotDataPack>();
            PlotData.Add(new ParentPlotDataPack("相移器电压(V)", "对比度数据", ChartDataType.X, Get1DChartDataSource("相移器电压(V)", "对比度数据"), false));//false有同名就不加，true加到结果里
            PlotData.Add(new ParentPlotDataPack("信号对比度[(sig-ref)/ref]", "对比度数据", ChartDataType.Y, Get1DChartDataSource("信号对比度[(sig-ref)/ref]", "对比度数据"), true));

            PlotData.Add(new ParentPlotDataPack("相移器电压(V)", "荧光数据", ChartDataType.X, Get1DChartDataSource("相移器电压(V)", "荧光数据"), false));
            PlotData.Add(new ParentPlotDataPack("平均光子数", "荧光数据", ChartDataType.Y, Get1DChartDataSource("平均光子数", "荧光数据"), true));
            PlotData.Add(new ParentPlotDataPack("信号光子数", "荧光数据", ChartDataType.Y, Get1DChartDataSource("信号光子数", "荧光数据"), true));
            return PlotData;
        }

        protected override List<KeyValuePair<string, Action>> AddPulseInteractiveButtons()
        {
            return new List<KeyValuePair<string, Action>>();
        }
    }
}
