﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Media;
using Controls.Charts;
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

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM.点实验.脉冲实验
{
    class T1 : PulseExpBase
    {
        public override string ODMRExperimentName { get; set; } = "驰豫时间测量(T1)";

        public override string ODMRExperimentGroupName { get; set; } = "点实验";

        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            new Param<int>("T1最小值(ns)",20,"T1min"),
            new Param<int>("T1最大值(ns)",100,"T1max"),
            new Param<int>("T1点数(ns)",20,"T1points"),
            new Param<int>("测量次数",1000,"LoopCount"),
            new Param<int>("序列循环次数",1000,"SeqLoopCount"),
            new Param<double>("微波频率(MHz)",2870,"RFFrequency"),
            new Param<double>("微波功率(dBm)",-20,"RFAmplitude"),
            new Param<int>("单点超时时间",10000,"TimeOut"),
        };

        public override List<ParamB> OutputParams { get; set; } = new List<ParamB>()
        {
        };

        public override List<KeyValuePair<DeviceTypes, Param<string>>> PulseExpDevices { get; set; } = new List<KeyValuePair<DeviceTypes, Param<string>>>();

        public override List<ChartData1D> D1ChartDatas { get; set; } = new List<ChartData1D>();
        public override List<ChartData2D> D2ChartDatas { get; set; } = new List<ChartData2D>();
        public override List<FittedData1D> D1FitDatas { get; set; } = new List<FittedData1D>();
        public override List<ODMRExpObject> SubExperiments { get; set; } = new List<ODMRExpObject>();

        public override bool IsAFMSubExperiment { get; protected set; } = true;

        protected override List<KeyValuePair<string, Action>> AddPulseInteractiveButtons()
        {
            return new List<KeyValuePair<string, Action>>();
        }

        public List<object> FirstScanEvent(object device, D1NumricScanRangeBase range, double locvalue, List<object> inputParams)
        {
            return ScanEvent(device, range, locvalue, inputParams);
        }

        public override bool PreConfirmProcedure()
        {
            if (MessageWindow.ShowMessageBox("提示", "是否要继续?此操作将清除原先的实验数据", MessageBoxButton.YesNo, owner: Window.GetWindow(ParentPage)) == MessageBoxResult.Yes)
            {
                return true;
            }
            return false;
        }

        private int CurrentLoop = 0;
        public List<object> ScanEvent(object device, D1NumricScanRangeBase range, double locvalue, List<object> inputParams)
        {
            //设置T1弛豫时间长度
            GlobalPulseParams.SetGlobalPulseLength("T1Step", (int)locvalue);

            PulsePhotonPack photonpack = DoPulseExp(GetInputParamValueByName("RFFrequency"), GetInputParamValueByName("RFAmplitude"), GetInputParamValueByName("SeqLoopCount"), 4, GetInputParamValueByName("TimeOut"));

            double signalcounts = photonpack.GetPhotonsAtIndex(0).Sum();
            double refcounts = photonpack.GetPhotonsAtIndex(1).Sum();

            int ind = range.GetNearestFormalIndex(locvalue);

            var contrfreq = Get1DChartDataSource("驰豫时间长度(ns)", "T1对比度数据");
            var signal = Get1DChartDataSource("对比度[(ref-sig)/ref]", "T1对比度数据");

            var florfreq = Get1DChartDataSource("驰豫时间长度(ns)", "T1荧光数据");
            var count = Get1DChartDataSource("参考光子数", "T1荧光数据");
            var sigcount = Get1DChartDataSource("信号光子数", "T1荧光数据");

            double signalcontrast = 0;
            try
            {
                signalcontrast = (signalcounts - refcounts) / refcounts;
            }
            catch (Exception)
            {
            }

            if (ind >= contrfreq.Count)
            {
                contrfreq.Add(locvalue);
                florfreq.Add(locvalue);
                signal.Add(signalcontrast);
                count.Add(refcounts);
                sigcount.Add(signalcounts);
            }
            else
            {
                signal[ind] = (signal[ind] * CurrentLoop + signalcontrast) / (CurrentLoop + 1);
                count[ind] = (count[ind] * CurrentLoop + refcounts) / (CurrentLoop + 1);
                sigcount[ind] = (sigcount[ind] * CurrentLoop + signalcounts) / (CurrentLoop + 1);
            }
            UpdatePlotChartFlow(true);
            return new List<object>();
        }

        public override void ODMRExpWithoutAFM()
        {
            int Loop = GetInputParamValueByName("LoopCount");
            double progressstep = 100 / Loop;
            for (int i = 0; i < Loop; i++)
            {
                CurrentLoop = i;
                Scan1DSession<object> Session = new Scan1DSession<object>();
                Session.FirstScanEvent = FirstScanEvent;
                Session.ScanEvent = ScanEvent;
                Session.ScanSource = null;
                Session.ProgressBarMethod = new Action<object, double>((obj, v) =>
                {
                    SetProgress(v);
                });
                Session.SetStateMethod = new Action<object, double>((obj, v) =>
                {
                    SetExpState("当前扫描轮数:" + i.ToString() + ",弛豫时间长度: " + Math.Round(v, 5).ToString());
                });

                D1NumricLinearScanRange range = new D1NumricLinearScanRange(GetInputParamValueByName("T1min"), GetInputParamValueByName("T1max"), GetInputParamValueByName("T1points"));

                Session.StateJudgeEvent = JudgeThreadEndOrResumeAction;
                Session.BeginScan(range, progressstep * i, progressstep * (i + 1));
            }
        }

        public override void PreExpEventWithoutAFM()
        {
            //打开微波
            RFSourceInfo RF = GetDeviceByName("RFSource") as RFSourceInfo;
            RF.Device.IsRFOutOpen = true;

            D1ChartDatas = new List<ChartData1D>()
            {
                new NumricChartData1D("驰豫时间长度(ns)","T1对比度数据",ChartDataType.X),
                new NumricChartData1D("对比度[(ref-sig)/ref]","T1对比度数据",ChartDataType.Y),

                new NumricChartData1D("驰豫时间长度(ns)","T1荧光数据",ChartDataType.X),
                new NumricChartData1D("参考光子数","T1荧光数据",ChartDataType.Y),
                new NumricChartData1D("信号光子数","T1荧光数据",ChartDataType.Y),
            };
            UpdatePlotChart();

            Show1DChartData("T1对比度数据", "驰豫时间长度(ns)", "对比度[(ref-sig)/ref]");
        }

        //T1拟合函数
        private double T1FitFunc(double x, double[] ps)
        {
            double a = ps[0];
            double tau = ps[1];
            double b = ps[2];
            return a * Math.Exp(-x / tau) + b;
        }

        public override void AfterExpEventWithoutAFM()
        {
            RFSourceInfo RF = GetDeviceByName("RFSource") as RFSourceInfo;
            RF.Device.IsRFOutOpen = false;
            //计算T1
            var xs = Get1DChartDataSource("驰豫时间长度(ns)", "T1对比度数据");
            var y1 = Get1DChartDataSource("对比度[(ref-sig)/ref]", "T1对比度数据");
            var tempdata = y1.Select(x => Math.Abs(x - (y1.Max() + y1.Min()) / 2)).ToList();
            double inittau = xs[tempdata.IndexOf(tempdata.Min())];
            double[] ps = CurveFitting.FitCurveWithFunc(xs, y1, new List<double>() { y1.Max() - y1.Min(), inittau, y1.Min() }, new List<double>() { 10, 10, 10 }, T1FitFunc, AlgorithmType.LevenbergMarquardt, 2000);

            //设置拟合曲线
            var ftxs = new D1NumricLinearScanRange(xs.Min(), xs.Max(), 500).ScanPoints;
            var fitys = ftxs.Select(x => T1FitFunc(x, ps)).ToList();
            D1FitDatas.Add(new FittedData1D("a*Exp(-x/t)+b", "x", new List<string>() { "a", "t", "b" }, ps.ToList(), "驰豫时间长度(ns)", "T1对比度数据", new NumricDataSeries("拟合曲线", ftxs, fitys) { LineColor = Colors.LightSkyBlue }));
            UpdatePlotChart();
            UpdatePlotChartFlow(true);
            Show1DFittedData("拟合曲线");
            OutputParams.Add(new Param<double>("T1拟合值(ns)", ps[1], "T1FitData"));
            //计算平均光子计数
            OutputParams.Add(new Param<double>("平均光子计数", Get1DChartDataSource("参考光子数", "T1荧光数据").Average(), "AverageCount"));
        }

        public override List<ParentPlotDataPack> GetD1PlotPacks()
        {
            List<ParentPlotDataPack> PlotData = new List<ParentPlotDataPack>();
            PlotData.Add(new ParentPlotDataPack("驰豫时间长度(ns)", "T1对比度数据", ChartDataType.X, Get1DChartDataSource("驰豫时间长度(ns)", "T1对比度数据"), false));
            PlotData.Add(new ParentPlotDataPack("对比度[(ref-sig)/ref]", "T1对比度数据", ChartDataType.Y, Get1DChartDataSource("驰豫信号对比度[sig]", "T1对比度数据"), true));

            PlotData.Add(new ParentPlotDataPack("驰豫时间长度(ns)", "T1荧光数据", ChartDataType.X, Get1DChartDataSource("驰豫时间长度(ns)", "T1荧光数据"), false));
            PlotData.Add(new ParentPlotDataPack("参考光子数", "T1荧光数据", ChartDataType.Y, Get1DChartDataSource("参考光子数", "T1荧光数据"), true));
            PlotData.Add(new ParentPlotDataPack("信号光子数", "T1荧光数据", ChartDataType.Y, Get1DChartDataSource("信号光子数", "T1荧光数据"), true));
            return PlotData;
        }

        protected override SequenceDataAssemble GetExperimentSequence()
        {
            return SequenceDataAssemble.ReadFromSequenceName("T1");
        }
    }
}
