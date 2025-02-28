﻿using ODMR_Lab.IO操作;
using ODMR_Lab.基本控件;
using ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM.点实验.CW谱扫描;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ODMR_Lab.实验部分.ODMR实验.实验方法.无AFM实验.单点.CW谱扫描
{
    internal class FourPointCW : CWBase
    {
        public override string ODMRExperimentName { get; set; } = "四频点CW";
        public override List<ParamB> InputParams { get; set; } = new List<ParamB>()
        {
            new Param<double>("频率点1",2850,"Frequency1"),
            new Param<double>("频率点2",2890,"Frequency2"),
            new Param<double>("频率点3",2850,"Frequency3"),
            new Param<double>("频率点4",2890,"Frequency4"),
            new Param<double>("微波功率(dBm)",-20,"RFPower"),
            new Param<int>("循环次数",1000,"LoopCount"),
            new Param<int>("单点扫描时间上限(ms)",0,"TimeOut"),
        };

        public override List<ParamB> OutputParams { get; set; } = new List<ParamB>();

        public override List<ParentPlotDataPack> GetD1PlotPacks()
        {
            return new List<ParentPlotDataPack>();
        }

        public override List<double> GetScanFrequences()
        {
            return new List<double>() { GetInputParamValueByName("Frequency1"), GetInputParamValueByName("Frequency2"),
            GetInputParamValueByName("Frequency3"), GetInputParamValueByName("Frequency4")};
        }

        protected override int GetLoopCount()
        {
            return GetInputParamValueByName("LoopCount");
        }

        protected override int GetPointTimeout()
        {
            return GetInputParamValueByName("TimeOut");
        }

        protected override double GetRFPower()
        {
            return GetInputParamValueByName("RFPower");
        }

        protected override void SetOutput()
        {
            List<double> contracts = GetContracts();
            List<double> counts = GetReferenceCounts();
            OutputParams.Add(new Param<double>("平均光子总数", counts.Average(), "Counts"));
            OutputParams.Add(new Param<double>("对比度1", contracts[0], "Contracts1"));
            OutputParams.Add(new Param<double>("对比度2", contracts[1], "Contracts2"));
            OutputParams.Add(new Param<double>("对比度3", contracts[2], "Contracts3"));
            OutputParams.Add(new Param<double>("对比度4", contracts[3], "Contracts4"));
        }
    }
}
