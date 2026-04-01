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
using ODMR_Lab.实验部分.ODMR实验.实验方法.其他.磁场调节.通用方法;
using ODMR_Lab.设备部分.位移台部分;
using ODMR_Lab.实验部分.ODMR实验.实验方法.其他;
using OpenCvSharp.XFeatures2D;
using System.Diagnostics.Contracts;
using ODMR_Lab.设备部分.电源;
//using System.Drawing;

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM.点实验.脉冲实验
{
    class PhaseShifterFine : PulseExpBase
    {
        public override bool Is1DScanExp { get; set; } = false;
        public override bool Is2DScanExp { get; set; } = true;

        public override string ODMRExperimentName { get; set; } = "相移器角度细调";

        public override string ODMRExperimentGroupName { get; set; } = "点实验";

        public override string Description { get; set; } = "使用的序列文件名称：PhaseShifter";

        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            new Param<double>("电压值(V)",0,"Voltage"),
            new Param<int>("测量次数",1000,"LoopCount"),        //外部的循环
            new Param<MultiScanType>("测量循环类型",MultiScanType.正向扫描,"ScanType"),
            new Param<int>("序列循环次数",1000,"SeqLoopCount"),   //板卡内的
            new Param<double>("微波频率(MHz)",2870,"RFFrequency"),
            new Param<double>("微波功率(dBm)",-20,"RFAmplitude"),
            new Param<int>("单点超时时间",10000,"TimeOut"),
        };
        public override List<ParamB> OutputParams { get; set; } = new List<ParamB>()
        {
        };
        public override List<KeyValuePair<DeviceTypes, Param<string>>> PulseExpDevices { get; set; } = new List<KeyValuePair<DeviceTypes, Param<string>>>() {
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
                if (D2ScanRange == null)
                {
                    MessageWindow.ShowTipWindow("扫描范围未设置", Window.GetWindow(ParentPage));
                    return false;
                }
                return true;
            }
            return false;
        }

        public List<object> ScanEvent(object arg1, object arg2, D2ScanRangeBase arg3, Point point, List<object> arg5)
        {
            //设置point
            GlobalPulseParams.SetGlobalPulseLength("HalfPiX", (int)point.X);
            GlobalPulseParams.SetGlobalPulseLength("HalfPiY", (int)point.Y);

            //大循环，记录总计数
            int loop = GetInputParamValueByName("LoopCount");
            double signalcount = 0;
            double refcount = 0;
            for (int i = 1; i <= loop; i++)
            {
                PulsePhotonPack pack = DoPulseExp("PhaseShifter", GetInputParamValueByName("RFFrequency"), GetInputParamValueByName("RFAmplitude"), GetInputParamValueByName("SeqLoopCount"), 4, GetInputParamValueByName("TimeOut"));
                signalcount += pack.GetPhotonsAtIndex(0).Sum();
                refcount += pack.GetPhotonsAtIndex(1).Sum();
            }
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
            int indx = D2ScanRange.GetNearestXIndex(point.X);
            int indy = D2ScanRange.GetNearestYIndex(point.Y);
            Get2DChartData("对比度", "扫描数据").Data.SetValue(indx, indy, signalcontrast);
            Get2DChartData("平均光子数", "扫描数据").Data.SetValue(indx, indy, refcount);
            Get2DChartData("信号光子数", "扫描数据").Data.SetValue(indx, indy, signalcount);
            UpdatePlotChart();
            UpdatePlotChartFlow(true);
            return new List<object>() { };
        }


        public override void ODMRExpWithoutAFM()
        {
            CustomScan2DSession<object, object> PointsScanSession = new CustomScan2DSession<object, object>();
            PointsScanSession.ProgressBarMethod = new Action<object, object, double>((d1, d2, v) => { SetProgress(v); });
            PointsScanSession.FirstScanEvent = ScanEvent;
            PointsScanSession.ScanEvent = ScanEvent;
            PointsScanSession.StartScanNewLineEvent = ScanEvent;
            PointsScanSession.EndScanNewLineEvent = ScanEvent;
            PointsScanSession.StateJudgeEvent = JudgeThreadEndOrResumeAction;
            PointsScanSession.SetStateMethod = new Action<object, object, Point>((d1, d2, p) =>
            {
                SetExpState("x脉冲时间:" + Math.Round(p.X, 5).ToString() + "  y脉冲时间:" + Math.Round(p.Y, 5).ToString());
            });
            PointsScanSession.BeginScan(D2ScanRange, 0, 100);
        }

        public override void PreExpEventWithoutAFM()
        {
            ////////////////rough
            SignalGeneratorChannelInfo RF = GetDeviceByName("RFSource") as SignalGeneratorChannelInfo;
            RF.Device.IsOutOpen = true;
            var power = GetDeviceByName("Power") as PowerChannelInfo;
            power.Channel.Voltage = GetInputParamValueByName("Voltage");

            //读取参数
            D1ChartDatas.Clear();
            D2ChartDatas.Clear();
            double xlo = D2ScanRange.XLo;
            double xhi = D2ScanRange.XHi;
            double ylo = D2ScanRange.YLo;
            double yhi = D2ScanRange.YHi;
            int xcount = D2ScanRange.XCount;
            int ycount = D2ScanRange.YCount;
            D2ChartDatas.Add(new ChartData2D(new FormattedDataSeries2D(xcount, xlo, xhi, ycount, ylo, yhi) { XName = "X脉冲时间", YName = "Y脉冲时间", ZName = "对比度" }) { GroupName = "扫描数据" });
            D2ChartDatas.Add(new ChartData2D(new FormattedDataSeries2D(xcount, xlo, xhi, ycount, ylo, yhi) { XName = "X脉冲时间", YName = "Y脉冲时间", ZName = "平均光子数" }) { GroupName = "扫描数据" });
            D2ChartDatas.Add(new ChartData2D(new FormattedDataSeries2D(xcount, xlo, xhi, ycount, ylo, yhi) { XName = "X脉冲时间", YName = "Y脉冲时间", ZName = "信号光子数" }) { GroupName = "扫描数据" });
            UpdatePlotChart();
            UpdatePlotChartFlow(true);
        }

        public override void AfterExpEventWithoutAFM()
        {
        }

        public override List<ParentPlotDataPack> GetD1PlotPacks()
        {
            return new List<ParentPlotDataPack>();
        }

        protected override List<KeyValuePair<string, Action>> AddPulseInteractiveButtons()
        {
            return new List<KeyValuePair<string, Action>>();
        }
    }
}
